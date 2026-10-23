#!/usr/bin/perl -w
use strict;

use Getopt::Long;

#############Befor  Start  , open the files ####################

sub usage
{
	print STDERR <<USAGE;
	2016-04-22       hewm\@genomics.cn

	Usage:		perl $0  -inPED in.ped -inMAP in.map  -outGenotype out.genotype

	Options
	-inPED        <s> :  Input the ped file path
	-inMAP        <s> :  Input the map file path
	-outGenotype  <s> :  Output the genotype path
	-maxMem       <i> :  Memory ceiling in MiB [half of MemAvailable]
	-block        <i> :  Variants per pass [derived from -maxMem]
	-verbose          :  Report the chosen mode on stderr
	-help             :  show this help

USAGE
}

my ($help,$inPED,$inMAP,$outGenotype,$optBlock,$optMaxMem,$optVerbose);

GetOptions(
	"help"=>\$help,
	"inPED:s"=>\$inPED,
	"inMAP:s"=>\$inMAP,
	"outGenotype:s"=>\$outGenotype,
	"block:i"=>\$optBlock,
	"maxMem:i"=>\$optMaxMem,
	"verbose"=>\$optVerbose,
);

if( defined($help) || !defined($outGenotype)   || !defined($inPED)   || !defined($inMAP)  )
{
	usage;
	exit ;
}


open (OA,">$outGenotype") || die "input file can't open $!";

my $Flag=6;
my %number2base=();
$number2base{"0"}="-";
$number2base{"1"}="A";
$number2base{"2"}="C";
$number2base{"3"}="G";
$number2base{"4"}="T";

# Allele code -> single base: numeric PED allele via %number2base, otherwise
# pass through unchanged (letters, "-", "N"). Undefined/empty becomes "-".
sub base_of
{
	my ($v)=@_;
	if(!defined($v) || $v eq "") { $v="-"; }
	return $number2base{$v} if exists $number2base{$v};
	return $v;
}

# Pair of allele codes -> one IUPAC token (native -InGenotype format: one
# token per sample). Same base -> homozygote single letter; two different
# bases -> heterozygote IUPAC (M/K/Y/R/W/S); missing pair -> "-"; mixed
# missing+base keeps the observed base.
my %het2iupac=("AC"=>"M","CA"=>"M","GT"=>"K","TG"=>"K","CT"=>"Y","TC"=>"Y",
               "AG"=>"R","GA"=>"R","AT"=>"W","TA"=>"W","CG"=>"S","GC"=>"S");
sub merge_alleles
{
	my ($x,$y)=@_;
	$x=base_of($x);
	$y=base_of($y);
	if($x eq "") { $x="-"; }
	if($y eq "") { $y="-"; }
	return $x if $x eq $y;
	return $y if $x eq "-";
	return $x if $y eq "-";
	my $pair="$x$y";
	return $het2iupac{$pair} if exists $het2iupac{$pair};
	return $pair;
}


# ---------------------------------------------------------------------------
# Memory. The previous version kept the entire PED in @Data as one scalar per
# field, i.e. nsamples x (6 + 2 nvariants) Perl scalars. On the paper's chr22
# panel that is 2,504 x 2,207,100 ~ 5.5e9 scalars, or 165-288 GiB, and it grows
# linearly with variants, so a whole-genome PED was impossible. The output only
# ever needs ONE merged IUPAC character per (sample, variant), which is one byte
# instead of ~40, so the data is streamed and the sample lines are discarded as
# soon as they have been folded into the output buffers.
#
# Buffers hold one string of up to nsamples characters per variant in the
# current block, so memory is nsamples x block bytes rather than nsamples x
# nvariants x 40. When the whole matrix does not fit the budget the PED is
# processed in blocks, re-reading it once per block: memory stays bounded, the
# cost is extra passes over the input.
# ---------------------------------------------------------------------------

# Variants, from the MAP file (one line per variant, no header).
my $nvar = 0;
{
	open(my $fm, "<", $inMAP) || die "input file can't open $!";
	$nvar++ while <$fm>;
	close $fm;
}
if ($nvar == 0) { die "$inMAP has no variants"; }

# Samples, from the PED file. Counted by scanning bytes rather than parsing, so
# this pass costs about one sequential read and no per-field work.
my $nsam = 0;
{
	open(my $fh, "<", $inPED) || die "input file can't open $!";
	binmode($fh);
	my ($buf, $n, $last);
	while (($n = read($fh, $buf, 1 << 20)) > 0) {
		$nsam += ($buf =~ tr/\n//);
		$last = substr($buf, -1, 1);
	}
	close $fh;
	# a final line without a trailing newline still counts as a sample
	$nsam++ if defined($last) && $last ne "\n";
}

# Memory budget in bytes.
#
# The default is chosen so that no tuning is needed: take half of what the
# kernel says is allocatable without swapping, which is enough to hold the whole
# matrix in the usual case and therefore to finish in a single pass - the
# fastest option, since blocking re-reads the PED once per block. A cap would
# only force pointless re-reads on machines that have the memory: the paper's
# chr22 panel needs 2 * 2504 * 1,103,547 bytes = 5.2 GiB, which an earlier
# 4 GiB cap would have pushed into block mode on a 502 GiB machine.
#
# -maxMem caps it explicitly for users who would rather bound memory than time;
# the block size is then derived, so there is nothing else to pass.
my $budget_mb = $optMaxMem;
if (!defined $budget_mb) {
    $budget_mb = 1024;                      # fallback when /proc is unreadable
    if (open(my $mi, "<", "/proc/meminfo")) {
        while (<$mi>) {
            if (/^MemAvailable:\s+(\d+) kB/) {
                $budget_mb = int($1 / 1024 / 2);   # half of what is free
                last;
            }
        }
        close $mi;
    }
    $budget_mb = 64 if $budget_mb < 64;      # never pick an absurd block
}
my $budget = $budget_mb * 1024 * 1024;

# The buffers are not the whole story: the interpreter, the line currently being
# read and the field array a fallback split produces all count towards RSS. A
# 200 MiB ceiling measured 260 MiB before this allowance was subtracted, i.e.
# the ceiling was not being honoured. Reserve a fixed amount and size the block
# from what is left, so -maxMem is a ceiling rather than a target.
my $OVERHEAD = 64 * 1024 * 1024;
# Never let the reservation eat the budget: with a 64 MiB ceiling and a 64 MiB
# reservation the usable part was zero, the block collapsed to one variant, and
# the PED was re-read once per variant -- thousands of passes, each re-splitting
# every line. Reserve at most half of the budget.
my $reserve = $OVERHEAD;
$reserve = int($budget / 2) if $reserve > int($budget / 2);
my $usable = $budget - $reserve;
$usable = 1 if $usable < 1;
if ($budget_mb < 32) {
	warn "warning: a $budget_mb MiB ceiling is below what the interpreter alone " .
	     "needs; the actual peak will exceed it. Use -maxMem 32 or more.\n";
}

my $block;
if (defined $optBlock) {
	die "-block must be at least 1\n" if $optBlock < 1;
	$block = $optBlock;
	$block = $nvar if $block > $nvar;
} elsif ($nsam > 0 && 2 * $nsam * $nvar <= $usable) {
	$block = $nvar;                       # one pass, whole matrix in buffers
} else {
	# two bytes per (sample, variant): the token plus its separator
	$block = int($usable / (2 * ($nsam > 0 ? $nsam : 1)));
	$block = 1 if $block < 1;
	$block = $nvar if $block > $nvar;
}
my $npass = int(($nvar + $block - 1) / $block);
if ($optVerbose) {
	printf STDERR "samples=%d variants=%d  block=%d  passes=%d  budget=%d MiB\n",
		$nsam, $nvar, $block, $npass, $budget_mb;
}
if ($npass > 1) {
	printf STDERR "note: processing %d variants in %d blocks to stay within the " .
		"%d MiB ceiling; the PED is read once per block. Raise -maxMem to trade " .
		"memory for speed.\n", $nvar, $npass, $budget_mb;
	if ($npass > 200) {
		printf STDERR "warning: %d passes will be slow; the cost grows with " .
			"(variants / block).\n", $npass;
	}
}

# Offset of the first genotype field: past the 6 mandatory PED columns. Returns
# -1 when the line has fewer than 6 fields, in which case the caller falls back
# to split.
sub geno_offset
{
	my ($line) = @_;
	my $len = length($line);
	my $off = 0;
	for (1 .. 6) {
		$off++ while $off < $len && substr($line, $off, 1) =~ /\s/;
		return -1 if $off >= $len;
		$off++ while $off < $len && substr($line, $off, 1) !~ /\s/;
	}
	$off++ while $off < $len && substr($line, $off, 1) =~ /\s/;
	return $off;
}

my $v0 = 0;                                # first variant of the current block
while ($v0 < $nvar) {
	my $v1 = $v0 + $block;
	$v1 = $nvar if $v1 > $nvar;
	my $nb = $v1 - $v0;

	# MAP lines for this block.
	my (@chr, @site);
	{
		open(my $fm, "<", $inMAP) || die "input file can't open $!";
		my $skip = $v0;
		<$fm> while $skip-- > 0;
		for (my $i = 0; $i < $nb; $i++) {
			my $line = <$fm>;
			last if !defined $line;
			chomp $line;
			my @inf = split ' ', $line;
			push @chr, $inf[0];
			push @site, abs($inf[-1]);
		}
		close $fm;
	}

	my @outbuf = ("") x $nb;
	my $uniform = -1;   # -1 unknown, 1 single-char fields, 0 irregular
	# Tokens are space separated, like the original join(" ", @geno). The
	# separator is added only from the second sample on, so the buffer holds
	# exactly the bytes the output line will contain.
	my $first_sam = 1;

	open(my $fp, "<", $inPED) || die "input file can't open $!";
	while (my $line = <$fp>) {
		chomp $line;

		# Two ways to reach the same allele pairs. Which one is cheaper
		# depends on how much of the line is needed:
		#
		#   single pass (nb == nvar)  split once and index. substr would
		#     allocate one scalar per field extracted, i.e. the same number
		#     of allocations as split, on top of split's own.
		#   block mode (nb < nvar)    extract only this block. Splitting
		#     first would parse every field of every line once per block,
		#     so M/K blocks would cost M/K * 2M extractions instead of 2M.
		#     A uniform PED allows plain substr; anything else falls back.
		if ($nb == $nvar) {
			my @inf = split ' ', $line;
			my $f = $Flag;
			for (my $i = 0; $i < $nb; $i++) {
				$outbuf[$i] .= ($first_sam ? "" : " ") .
					merge_alleles($inf[$f], $inf[$f + 1]);
				$f += 2;
			}
		} else {
			my $fast = 0;
			my $off = geno_offset($line);
			if ($off >= 0) {
				my $rem = substr($line, $off);
				# A PED has one line per sample with the same layout, so
				# the structural check runs on the first line only; later
				# lines are accepted on length. A line of a different
				# length falls back, which is what ragged input needs.
				if (length($rem) == 4 * $nvar - 1) {
					if ($uniform < 0) {
						$uniform = ($rem =~ /^\S(?: \S)*$/) ? 1 : 0;
					}
					$fast = 1 if $uniform;
					$fast = 1 if !$uniform && $rem =~ /^\S(?: \S)*$/;
				}
				if ($fast) {
					for (my $i = 0; $i < $nb; $i++) {
						my $p = 4 * ($v0 + $i);
						$outbuf[$i] .= ($first_sam ? "" : " ") .
							merge_alleles(substr($rem, $p, 1),
							              substr($rem, $p + 2, 1));
					}
				}
			}
			if (!$fast) {
				my @inf = split ' ', $line;
				my $f = $Flag + 2 * $v0;
				for (my $i = 0; $i < $nb; $i++) {
					$outbuf[$i] .= ($first_sam ? "" : " ") .
						merge_alleles($inf[$f], $inf[$f + 1]);
					$f += 2;
				}
			}
		}
		$first_sam = 0;
	}
	close $fp;

	for (my $i = 0; $i < $nb; $i++) {
		my $c = defined $chr[$i] ? $chr[$i] : "";
		my $s = defined $site[$i] ? $site[$i] : 0;
		print OA "$c\t$s\t$outbuf[$i]\n";
	}
	$v0 = $v1;
}
close OA ;

######################swiming in the sky and flying in the sea #############################

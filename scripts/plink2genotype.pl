#!/usr/bin/perl-w
use strict;

use Data::Dumper;
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

	-help             :  show this help

USAGE
}

my ($help,$inPED,$inMAP,$outGenotype);


GetOptions(
	"help"=>\$help,
	"inPED:s"=>\$inPED,
	"inMAP:s"=>\$inMAP,
	"outGenotype:s"=>\$outGenotype,
);

if( defined($help) || !defined($outGenotype)   || !defined($inPED)   || !defined($inMAP)  )
{
	usage;
	exit ;
}


open (IA,"$inPED") || die "input file can't open $!";
open (IB,"$inMAP") || die "input file can't open $!";
open (OA,">$outGenotype") || die "input file can't open $!";

my @Data;

my $sampleNum=0;
while(<IA>)
{
	chomp ;
	my @inf=split ;
	push @Data,\@inf;
	$sampleNum++;
}
close IA;
$sampleNum--;
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


while(<IB>)
{
	chomp ;
	my @inf=split ;
	my $chr=$inf[0];
	my $site=abs($inf[-1]);
	my $SedFlag=$Flag+1;
	my @geno;
	for my $samID (0..$sampleNum)
	{
		my $A=$Data[$samID]->[$Flag];
		my $B=$Data[$samID]->[$SedFlag];
		push @geno, merge_alleles($A,$B);
	}

	$Flag++;
	$Flag++;
	print OA "$chr\t$site\t".join(" ", @geno)."\n";
}
close IB ;
close OA ;

######################swiming in the sky and flying in the sea #############################

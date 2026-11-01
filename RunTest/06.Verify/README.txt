06.Verify -- the regression gate for the binary you just built
===============================================================

Run this after ./make.sh to check that PopLDdecay2 still reproduces its expected
output. It needs no network, no reference data and no PopLDdecay 3.45, so this
is also what continuous integration runs.

  bash step1_checkOutput.sh          ### 108 cases, prints PASS/FAIL, exits non-zero on failure

What it does: regenerates the inputs with the seeded generators in this
directory, runs bin/PopLDdecay2 on every command line in expected.tsv, and
compares the SHA-256 of the output (after gunzip for .stat.gz) against the
frozen digest. It then runs the behaviour cases, which pin the exit status and
the message of a run that cannot succeed (an invalid number, a missing
argument, an unwritable output directory) instead of a digest, and the
documentation cases, which keep README.md and README.cn.md from drifting apart
(the same option names, working links, the case count).

  expected.tsv        ### 88 frozen cases: case_id, kind, argv, suffix, sha256
  gen_inputs.sh       ### builds the synthetic inputs, shared by both scripts below
  gen_sim_vcf.py      ### seeded VCF generator
  gen_sim_genotype.py ### VCF -> -InGenotype
  gen_sim_ped.py      ### PED/MAP generator, the input of plink2genotype.pl
  gen_sim_plink.py    ### PLINK .bed/.bim/.fam generator
  golden_diff.sh      ### byte-for-byte comparison against PopLDdecay 3.45
  make_fixtures.sh    ### regenerate expected.tsv after an intentional change

The 83 digests of the previous release are unchanged: the five new ones (a
200-character chromosome name, whose rows used to overflow the buffers of the
.LD.gz and .ehh.gz writers) were added on top.

golden_diff.sh is the only script here that needs PopLDdecay 3.45: put it in the
package root as PopLDdecay-3.45/bin/PopLDdecay and run

  bash golden_diff.sh

It uses the same case matrix as expected.tsv but compares against the real old
binary instead of frozen digests, so it needs that binary to say anything. It
also runs the cases that need the public 1000 Genomes chr22
(RunTest/00.GetData/chr22/chr22.vcf.gz) and a real plink2; they are skipped,
each with a line saying what is missing, when those are not present.

Regenerate expected.tsv only when an output change is intentional: commit the new
file together with the code change and say why in the commit message.

Note:
1. bgzip and tabix from samtools must be in PATH; step1_checkOutput.sh refuses to start
   without them.
2. step1_checkOutput.sh builds its inputs in a temporary directory and removes it on exit.
   Set CI_WORK=/path to keep them.
3. This directory was called 06.SelfTest before; the scripts and the paths in
   the documentation were renamed with it.


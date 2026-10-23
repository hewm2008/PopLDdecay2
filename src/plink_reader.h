#pragma once

#include <string>
#include <vector>

#include "genotype_reader.h"  // GenoSite / GenoBuf
#include "subpop.h"
#include "types.h"

namespace pld2 {

// Reads a PLINK binary dataset (.bed + .bim + .fam) into the same unpacked
// per-contig buffers the genotype reader produces, so `merge_geno_buffers`
// finishes the job exactly as it does for -InGenotype.
//
// `path` is either the shared prefix ("data/chr22") or the .bed file itself
// ("data/chr22.bed"); .bim/.fam are derived from it. Samples are the .fam IIDs
// (column 2), matching what `plink --recode vcf` puts in the VCF header.
//
// Site filters and encoding replicate the VCF reader, not the genotype reader,
// so that for one dataset `-InPlink x.bed` and `-InVCF x.vcf` are byte
// identical. Variants whose .bim alleles are not single A/C/G/T bases are
// dropped by the Indel filter; a1 == a2 is dropped as non-bi-allelic.
//
// -SubPop, -L, -OutFilterSNP and the parallel chunking behave as they do for
// the other readers; see plink_reader.cpp for the details. Returns 1 on
// success, 0 on error (with a message).
int read_plink(const std::string& path, const Options& opt,
               const std::string& prefix,
               const std::vector<SubGroupSpec>& specs, ReadStats& stats,
               std::vector<std::vector<GenoBuf>>& out);

}  // namespace pld2
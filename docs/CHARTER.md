# Charter

just-vcf turns genotype data that someone else already produced into VCF: consumer array exports,
freeform genotype tables, and VCF/gVCF from any caller. This document holds what every release
upholds. It outranks any plan.

## What it is

- A converter and normalizer. It reads a genotype file, works out what kind of file it is and what it
  can honestly say about each position, and writes that out as VCF.
- Built on polars-bio and Polars, lazy wherever the input is large.

## What it will never do

- Call variants from reads. BAM, CRAM and FASTQ are out of scope.
- Impute genotypes. Filling a position from population frequencies or haplotype panels is a different
  tool's job, and output from such a tool is an input here like any other.
- Present its output as clinical-grade.

## Invariants

1. **Genotype evidence is graded, never binary.** A position is not "called" or "not called". Its
   evidence runs from none (unknown) to a confident genotype, on a scale like a PHRED quality. The
   representation is open. It is probably a vector of genotype probabilities or likelihoods per
   position. Whatever it becomes, no public shape reduces it to a boolean. Inputs that disagree with
   each other (several callers, imputed sets) are expected. They are not an error.
2. **Unknown is a value.** Anything answering a question (which build, which format, which callset
   type, which population) can answer "unknown", and says so rather than guessing a default.
3. **REF comes from a reference.** REF is read from a reference sequence or a reference panel, or it is
   left unresolved. It is never inferred from the sample's own alleles.
4. **Absence is evidence only where the input type supports it.** A site missing from a whole-genome
   callset is weak evidence for reference. A site missing from an array export is no evidence at all.
   What an absence means is decided by the detected input type and carried as graded evidence (1).
5. **Every output position says where it came from**: the input file, the reference, or a
   restoration rule.
6. **Output is deterministic.** The same input, reference and options give byte-identical output.
7. **Coordinates at every public boundary are 1-based VCF positions.**

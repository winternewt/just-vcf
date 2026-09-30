# Charter

just-vcf holds the results of variant calling for a person as Arrow data, read through Polars and
polars-bio. VCF and gVCF in the name are a tribute. The files are what comes in, but the product is
the dataframe. This document holds what every release upholds. It outranks any plan.

## What it is

- **Minimal shape: many formats in, one schema out.** VCF, gVCF, consumer array exports and freeform
  genotype tables are read into three normalized, lazy frames: `facets` (what each input is and how it
  was made), `calls` (one row per site) and `blocks` (reference and no-call intervals). A gVCF does not
  fit in a calls frame alone. Writing files back out, VCF included, is an export at the edge.
- **Extended shape: every representation of a genome, side by side and quantified.** One person can
  have calls from different assemblies (GRCh37, GRCh38, T2T-CHM13), different aligners and liftover
  routes, different callers (DRAGEN, GATK, DeepVariant, Clair3, Dorado), and different phasings (HiFi
  reads, trio, imputation). Coverage is a facet too, as evidence of what was callable: summaries from
  samtools or mosdepth, or depth computed from alignments by what polars-bio already provides. Each of these is a **facet**. just-vcf keeps them all, says how each one
  was produced, measures where they agree and disagree, and links them in a graph. The graph's shape is
  open. It starts as links between facets and alleles. A sequence or variation graph over the alleles
  themselves is neither planned nor ruled out, and nothing here may close that path.
- **Arrow-native, in the CompoSeq way.** Like polars-bio and oxbow, it works on columnar data with no
  custom in-memory object model, streams wherever the input is large, and does not re-parse what an
  existing Arrow reader already reads.
- **GA4GH from the start.** Allele identity is VRS. Sequence and assembly identity is refget and refget
  Sequence Collections. Phasing is VRS CisPhasedBlock. Pedigrees, when they come, are the GA4GH
  Pedigree standard. Cross-assembly links use the Cat-VRS `liftover_to` relation.

## What it will never do

- Implement alignment processing. just-vcf writes no BAM, CRAM or FASTQ handling of its own. Where the
  substrate already provides something over alignments (polars-bio's depth), it is used as is.
  Aligners, callers, phasers and imputation servers make the facets, and just-vcf reads what they
  produced.
- Present its output as clinical-grade.
- Emit a GA4GH-shaped value that does not conform. An identifier that looks interoperable and is not is
  worse than no identifier.

## Invariants

1. **The frame is the product.** Every public operation takes and returns Arrow/Polars data. A file
   format is an input or an export, never the internal model.
2. **GA4GH first, our own shape only where GA4GH has none.** Where a GA4GH standard defines a model, an
   identifier or a term, just-vcf uses it and its vocabulary verbatim. Where it defines none, the shape
   is just-vcf's own and is documented as such, never dressed in GA4GH names.
3. **An assembly is identified by its content.** A facet's reference is identified by its refget
   Sequence Collection digest. "GRCh38", "hg38" and "T2T-CHM13" are labels on top of that. hg38 with
   alt contigs, the no-alt analysis set and T2T are different collections, and the digest is what tells
   them apart.
4. **An allele is identified by VRS.** A call's identity is the VRS Allele id of its fully justified
   allele on a refget-addressed sequence. That id is the join key between facets on the same assembly.
   When an id cannot be minted conformantly (an indel that cannot be justified, an assembly with no
   refget table yet), it is withheld, not approximated, and the row stays usable without it.
5. **Nothing is overwritten. Facets accumulate.** Importing a new callset adds a facet; it does not
   replace one. A derived facet keeps a link to what it came from: `liftover_to` for a liftover, and a
   documented just-vcf relation for anything GA4GH has no term for (a merge, a re-call, a consensus).
   When a more direct facet arrives, such as calls from reads aligned to GRCh37 alongside calls lifted
   from GRCh38, the direct one is **preferred**. The derived one stays.
6. **Provenance is part of every row.** Each call can be traced to its facet, and each facet records its
   method: assay, reference, aligner, caller and version, phasing method, and any post-processing.
7. **Genotype evidence is graded, never binary.** A position is not "called" or "not called". Its
   evidence runs from none (unknown) to a confident genotype, on a scale like a PHRED quality.
   VRS 2.x has no genotype model, so this representation is just-vcf's own (invariant 2). It is
   probably a vector of genotype probabilities or likelihoods per allele set. No public shape reduces
   it to a boolean.
8. **Disagreement is data.** Concordance and discordance between facets are measured and kept, never
   resolved away by silently picking a winner. Systematic errors are attributed to the method that
   makes them, so a facet's reliability is stated by method and by region, not assumed.
9. **Unknown is a value.** Anything answering a question (which assembly, which format, which caller,
   which population) can answer "unknown", and says so rather than guessing a default.
10. **REF comes from the reference.** REF is read from the facet's own assembly sequence, addressed by
    its refget accession, or left unresolved. It is never inferred from the sample's alleles.
11. **What an absent position means depends on the facet.** The facet's type (whole-genome, exome,
    array, gVCF blocks, imputed) decides it, and the result is carried as graded evidence (7).
12. **Output is deterministic.** The same inputs, references and options give identical data.
13. **Two coordinate systems, never mixed.** Plain position columns are 1-based VCF positions. VRS
    objects keep the standard's 0-based interbase coordinates. They never share a column name, and the
    conversion between them happens in one place and is tested.

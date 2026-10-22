# General

Stratification group: XY

These stratifications represent the sex chromosomes.

For contact information, sharing and access information, data usage policy,
recommended citations, general workflow information, quality control, and other
high-level information, see the top-level readme one level up.

# File Descriptions

## All Autosomes

This represents all autosomal regions (minus gaps).

`GRCh37_AllAutosomes.bed.gz`

## Sex Chromosome Files

X and Y chromosomes are broken out into separate files denoted with "chrX" or
"chrY".

### Pseudo-autosomal Regions (PAR)

These are regions at both ends that are homologous between X and Y.

`GRCh37_chrX_PAR.bed.gz`
`GRCh37_chrY_PAR.bed.gz`


The complements of these (nonPAR) are given in:

`GRCh37_chrX_nonPAR.bed.gz`
`GRCh37_chrY_nonPAR.bed.gz`


### Ampliconic

These are regions on both X and Y that tend to have many copies.

`GRCh37_chrX_ampliconic.bed.gz`
`GRCh37_chrY_ampliconic.bed.gz`


### X-transposed

These are regions that were transposed from X to Y and thus are highly
homologous between the two (but also outside the PAR).

`GRCh37_chrX_XTR.bed.gz`
`GRCh37_chrY_XTR.bed.gz`


# Provenance and Data Processing

## All autosomes

The autosomes bed file was generated using the entire length of chromosomes 1-22
and subtracting gaps.

## PAR

Each of these was used as-is and also complemented with `complementBed` to
produce the non-PAR beds.

### X

The PAR coordinates were specified manually in the pipeline configuration.

### Y

The PAR coordinates were specified manually in the pipeline configuration.

## Features

### X

The X chromosome feature bed was downloaded from
https://raw.githubusercontent.com/SexChrLab/SexChrCoordinates/4ad1faf2b822db6fb89498a750d0b5336703fde3/hg19/chrX_genomic_features_hg19.bed.
The md5 hash of the uncompressed file was bbe9bf0948e1c87d9aaf22aca432b52e.

To construct the bed file, columns 1, 2, and 3 were selected as the 'chrom',
'start', and 'end' columns respectively.

Specific features were isolated by filtering column 4 for 'XTR' and 'Ampliconic'

### Y

The Y chromosome feature bed was downloaded from
https://raw.githubusercontent.com/SexChrLab/SexChrCoordinates/4ad1faf2b822db6fb89498a750d0b5336703fde3/hg19/chrY_genomic_features_hg19.bed.
The md5 hash of the uncompressed file was d22c3c3c73c9671fc1bcf5914af867ce.

To construct the bed file, columns 1, 2, and 3 were selected as the 'chrom',
'start', and 'end' columns respectively.

Specific features were isolated by filtering column 4 for 'XTR' and 'Ampliconic'

## Software

- bedtools: version 2.30.0
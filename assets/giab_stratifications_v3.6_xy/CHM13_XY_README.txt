# General

Stratification group: XY

These stratifications represent the sex chromosomes.

For contact information, sharing and access information, data usage policy,
recommended citations, general workflow information, quality control, and other
high-level information, see the top-level readme one level up.

# File Descriptions

## All Autosomes

This represents all autosomal regions (minus gaps).

`CHM13_AllAutosomes.bed.gz`

## Sex Chromosome Files

X and Y chromosomes are broken out into separate files denoted with "chrX" or
"chrY".

### Pseudo-autosomal Regions (PAR)

These are regions at both ends that are homologous between X and Y.

`CHM13_chrX_PAR.bed.gz`
`CHM13_chrY_PAR.bed.gz`


The complements of these (nonPAR) are given in:

`CHM13_chrX_nonPAR.bed.gz`
`CHM13_chrY_nonPAR.bed.gz`


### X-transposed

These are regions that were transposed from X to Y and thus are highly
homologous between the two (but also outside the PAR).

`CHM13_chrX_XTR.bed.gz`
`CHM13_chrY_XTR.bed.gz`


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
https://raw.githubusercontent.com/SexChrLab/SexChrCoordinates/4ad1faf2b822db6fb89498a750d0b5336703fde3/T2T/chrX_genomic_features_CHM13.bed.
The md5 hash of the uncompressed file was 4a49d47fff4a3fac778ec66aff5c19b1.

To construct the bed file, columns 1, 2, and 3 were selected as the 'chrom',
'start', and 'end' columns respectively.

Specific features were isolated by filtering column 4 for 'XTR'

### Y

The Y chromosome feature bed was downloaded from
https://raw.githubusercontent.com/SexChrLab/SexChrCoordinates/4ad1faf2b822db6fb89498a750d0b5336703fde3/T2T/chrY_genomic_features_HG002.bed.
The md5 hash of the uncompressed file was 5b78ffe218c11c18b6a0eaeb0c7fdf76.

To construct the bed file, columns 1, 2, and 3 were selected as the 'chrom',
'start', and 'end' columns respectively.

Specific features were isolated by filtering column 4 for 'XTR'

## Software

- bedtools: version 2.30.0
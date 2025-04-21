# Genotype PRDM9 from long reads
Script to genotype PRDM9 Zinc-finger array from reads that span the entire locus. Adapted from Alleva et al. 2021. 

## Requirements:
* Singularity 3.6.4+ or docker 27.5.1+

## Notes:
Can use long reads from any current sequencing technology. See Alleva et al. 2021 for details.

### Installation: 
Download genotypePRDM9FromLongReads to a location on your path.

### Run :
genotypePRDM9FromLongReads 

#### Usage
```
usage: genotypePRDM9 -f <<fasta file>> -n <<name>> [-v] [-h]
-------------------------------------------------------------------
Genotype PRDM9 from long sequencing reads
(See Alleva et al. 2021)
-------------------------------------------------------------------
-f/--fa      : FASTA file containing long reads from
               the PRDM9 ZF-Array (Required)

-n/--name    : Sample name (Required)

-a/--alleles : File of PRDM9 alleles
               All PRDM9 alleles used in Alleva et al. 2021 are
               queried by default. This option is useful if you
               wish to include new alleles. File format must be as
               found in File S3 of Alleva et al. 2021

-z/--zfs     : File of PRDM9 ZFs
               All PRDM9 ZFs used in Alleva et al. 2021 are
               queried by default. This option is useful if you
               wish to include new ZFs. File format must be as
               found in File S2 of Alleva et al. 2021

-t           : Run internal test

-v           : verbose output files

-d/--docker  : Use Docker instead of Singularity

-h           : show this help
-------------------------------------------------------------------
Outputs:
PRDM9 haplotypes table      : <<name>>.PRDM9alleles.text
                              One line per allele

ZF arrays from long reads   : <<name>>.PrZFA.##ZFs.fa
                              One file for each size ZF-array
-------------------------------------------------------------------
**** ARGUMENT ERROR ****
Sample name required: -n/--name <<name>>
```
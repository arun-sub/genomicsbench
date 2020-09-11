Three kernels from PLINK are taken.

### Datasets needed:

PLINK1.9: Download chr1 dataset from [here](https://www.cog-genomics.org/plink/2.0/resources#1kg_phase3). We need the .pgen, .pvar and pedigree-corrected .psam files.
PLINK2: Download chr1 dataset from [here](ftp://parrot.genomics.cn/gigadb/pub/10.5524/100001_101000/100116/1kg_phase1_chr1.tar.gz).

### Command lines
* --make-grm-bin (PLINK1.9 and PLINK2)

1.9: plink --bfile <1kg\_phase1\_chr1> --make-grm-bin --out <outfile prefix>
2: plink2 --maf 0.01 --pgen <> --pvar <> --psam <> --make-grm-bin --out <outfile prefix>

* --genome (PLINK1.9)
plink --noweb --bfile <> --genome --out <outfile prefix>

* --epistasis (PLINK1.9)
plink --noweb --bfile UC\_N6000\_M7000\_P1\_C5/plink --epistasis --out <outfile prefix>

# PanTree Contributor Notes

PanTree 0.0.1 is an analysis tool for VCFs generated from minigraph cactus's panegenome tool for assessing variant
occupancy and phylogenetic relationship among species and breed-level pangenomes.

General: supplied files are  VCFs and (optionally) reference genome FASTAs
- VCFs represent whole value genotypes (0, 1, 2, etc.).
- VCFs have been pre-filtered to remove any variants with Ns in genotypes
- VCFs have been pre-filtered for preferred missing genotype calls and variant states 

Core functions:

- Read in VCF
- Chunk genome/chromosomes into windows (by basepair windows, variable variant count; and by variant count, variable
window size)
- Calculate the following across windows:
  - called genotype density
- Build sequence alignments based on alleles
  - SNPs and MNPs treated as prealigned sites
  - Indels must be aligned via MAFFT (Biopython)
  - Multi-state variants must be aligned via MAFFT (Biopython)
- Tabulate information for window
  - Chromosome
  - Window start bp
  - Window stop bp
  - Total variant sites
  - Total SNPs
  - Total MNPs
  - Total INDELs
  - Total basepairs
  - Neighbor-joining tree from window alignment
  - For each sample
    - Total sites with called variant
    - Total variants SNPs
    - Total variants MNPs
    - Total variants INDELs
    - Total basepairs
    - Total unique alleles (e.g., 1 versus all others 0; 0 permissible if Reference is not considered, determined
    in CLI)
  - For each sample to sample pair (samples of interest defined in CLI)
    - Patristic (tip to tip) distance based onb NJ tree
    - Jaccard distances (only calculated if Reference is not considered)
    - Manhattan distance normalized by valid sites

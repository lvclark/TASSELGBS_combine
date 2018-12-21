# TASSELGBS_combine
Pipeline to combine SNP datasets from multiple TASSEL-GBS databases

This is a pipeline for a collaboration between the Sacks and Lipka labs, where we are 
experimenting with performing genomic prediction across different diversity panels and
mapping populations.  SNP mining has been performed on each of these populations 
independently using the TASSEL-GBSv2 pipeline.  For a given set of populations, we want
to identify RAD-seq tag locations (i.e. restriction cut site * strand) that were present 
in all populations.  Then we want to filter markers to retain in the combined dataset,
where our requirements for missing data rate and minor allele frequency might vary from
population to population.  The final output will be an accession x allele genotype matrix, 
with numbers ranging from 0 to 1 indicating allele copy number divided by the ploidy of 
each individual.

## Pipeline outline

1. In TASSEL-GBSv2, for each population, run `GBSSeqToTagDBPlugin` and `TagExportToFastqPlugin`.
For the latter, it is advisable to set `-c` to something higher than the default.
2. Use Bowtie2 or BWA to generate one SAM file per population, using the FASTQ files created
by the above step.
3. A script using [TagDigger](https://github.com/lvclark/TagDigger) functions identifies common 
tag locations across SAM files, organizes tags into markers (by location), and exports tag 
sequences in TagDigger's database format.
4. TagDigger is used for obtaining read depth of each selected tag in each individual across
all populations.
5. The read depth matrix is evaluated for missing data rate and minor allele frequency (or 
number of individuals with minor allele) per population, which is used for filtering markers.
6. The read depth matrix is split by population * ploidy, and genotypes are called using 
[polyRAD](https://github.com/lvclark/polyRAD).
7. Genotypes are exported on a scale of 0 to 1 using `GetWeightedMeanGenotypes` on each
population * ploidy, combined into one matrix with `rbind`, then exported to CSV.

# Merged phenotypes

Do we want to use LS means to get phenotypic values for the entire combined set all at once
with population (i.e. experiment) as an effect?

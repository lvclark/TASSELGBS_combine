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

1. In TASSEL-GBSv2, for each population, run `GBSSeqToTagDBPlugin`,
`GetTagTaxaDistFromDBPlugin`, and `TagExportToFastqPlugin`.
For the latter, it is advisable to set `-c` to something higher than the default.
2. Use Bowtie2 or BWA to generate one SAM file per population, using the FASTQ files created
by the above step.
3. Edit `params.txt` to indicate where the files are located, what SNP filtering
parameters you would like to use, and where files should be output.
4. Run `common_tags.py` using Python 3. Using [TagDigger](https://github.com/lvclark/TagDigger)
functions, this script identifies common tag locations across SAM files and
organizes tags into markers (by location).  The read depth of each selected tag
in each individual across all populations is retrieved from the TagTaxaDist table
exported from the TASSEL-GBS database in step 1.  The read depth matrix is
evaluated for missing data rate and number of individuals with minor allele per
population, which is used for filtering markers. For this final set of markers,
tag sequences are exported in TagDigger's database format.
5. Edit `call_genotypes.R` and run it using R 3.5. The read depth matrix is split by 
population * ploidy, and genotypes are called using
[polyRAD](https://github.com/lvclark/polyRAD). Genotypes are exported on a scale
of 0 to 1 using `GetWeightedMeanGenotypes` on each
population * ploidy, combined into one matrix with `rbind`, then exported to CSV.

# Merged phenotypes

Do we want to use LS means to get phenotypic values for the entire combined set all at once
with population (i.e. experiment) as an effect?

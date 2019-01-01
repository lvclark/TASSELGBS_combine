# Take the read counts exported by common_tags.py and call genotypes using polyRAD

library(polyRAD)
library(qqman)

paramsFile <- "params.txt" ## Edit this file name if necessary

# import parameters, including file and directory names
paramsLines <- readLines(paramsFile)

getParams <- function(linestart, paramsl = paramsLines){
  val <- grep(linestart, paramsl, value = TRUE)
  val <- trimws(sub(linestart, "", val))
  return(val)
}

tagdb_file <- getParams("Tag database output file:")
outdir <- getParams("Output directory:")
counts_files <- file.path(outdir , paste(getParams("Population name:"), 
                                         "_counts.csv", sep = ""))

## Genotype calling; edit as necessary for different populations

# Call genotypes in the test population first; this is the one where we care the
# most about alleles being variable.

# 09F2 population ####
RD_09F2 <- readTagDigger(grep("09F2", counts_files, value = TRUE),
                         dbfile = tagdb_file, possiblePloidies = list(2),
                         dbChrCol = "Chromosome", dbPosCol = "Position")
RD_09F2 <- SubsetByTaxon(RD_09F2, 
                         GetTaxa(RD_09F2)[!GetTaxa(RD_09F2) %in% c("DDH2O", "H2O", "blank")])
RD_09F2 <- SetDonorParent(RD_09F2, "UI10-00009")
RD_09F2 <- SetRecurrentParent(RD_09F2, "UI10-00014")

RD_09F2 <- AddPCA(RD_09F2)
plot(RD_09F2) # all seems to be one population

RD_09F2_prelim <- PipelineMapping2Parents(RD_09F2, n.gen.intermating = 1,
                                   useLinkage = FALSE)

odtest <- TestOverdispersion(RD_09F2_prelim, to_test = 8:10)
qq(odtest[[1]])
qq(odtest[[2]])

RD_09F2 <- PipelineMapping2Parents(RD_09F2, n.gen.intermating = 1,
                                   useLinkage = TRUE, overdispersion = 8)

mat_09F2 <- GetWeightedMeanGenotypes(RD_09F2, omit1allelePerLocus = FALSE)

# preview results
unname(round(mat_09F2[1:10,1:10], 2))
sum(is.na(mat_09F2))/nrow(mat_09F2)
freq <- colMeans(mat_09F2, na.rm = TRUE)
hist(freq, breaks = 50)
sum(freq > 0.01 & freq < 0.99, na.rm = TRUE)

# Subset to things that are actually segregating in the pop, allowing for a 
# lot of distortion.
mat_09F2 <- mat_09F2[, which(freq > 0.01 & freq < 0.99)]
dim(mat_09F2) # 293 ind x 4766 alleles

# visualize
head(colnames(mat_09F2))
chr01 <- grep("Chr01", colnames(mat_09F2))
image(mat_09F2[,chr01[seq(2, length(chr01), by = 2)]])

#save(mat_09F2, file = file.path(outdir, "mat_09F2.RData"))
#write.csv(mat_09F2, file = file.path(outdir, "genotype_mat_09F2.csv"))

load(file.path(outdir, "mat_09F2.RData"))

# alleles per marker?
mrkr <- sub("_[ACGTN]*$", "", colnames(mat_09F2))
table(table(mrkr)) # mostly 2, some 3 or 4, some unexpected

hist(mat_09F2["UI10-00009",])
hist(mat_09F2["UI10-00014",]) # both parents have a lot of heterozygous genotypes

# Msa diploids ####

# Msa tetraploids ####

# Msi population ####

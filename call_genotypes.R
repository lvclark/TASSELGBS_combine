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

# 09F2 population with PipelineMapping2Parents; not used ####
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

# Alternative 09F2 under HWE model ####
# I think this may better deal with segregation distortion
RD_09F2 <- readTagDigger(grep("09F2_counts", counts_files, value = TRUE),
                         dbfile = tagdb_file, possiblePloidies = list(2),
                         dbChrCol = "Chromosome", dbPosCol = "Position")
RD_09F2 <- SubsetByTaxon(RD_09F2, 
                         GetTaxa(RD_09F2)[!GetTaxa(RD_09F2) %in% c("DDH2O", "H2O", "blank")])

RD_09F2_prelim <- IterateHWE(RD_09F2)

odtest <- TestOverdispersion(RD_09F2_prelim, to_test = 6:10)
# tiff(file.path(outdir, "09F2overdispersion6.tiff"), compression = "lzw")
# qq(odtest[["6"]])
# dev.off()
# tiff(file.path(outdir, "09F2overdispersion7.tiff"), compression = "lzw")
# qq(odtest[["7"]])
# dev.off()
# tiff(file.path(outdir, "09F2overdispersion8.tiff"), compression = "lzw")
# qq(odtest[["8"]])
# dev.off()
# tiff(file.path(outdir, "09F2overdispersion9.tiff"), compression = "lzw")
# qq(odtest[["9"]])
# dev.off()
# tiff(file.path(outdir, "09F2overdispersion10.tiff"), compression = "lzw")
# qq(odtest[["10"]])
# dev.off()

rm(RD_09F2_prelim)

RD_09F2 <- IterateHWE_LD(RD_09F2, LDdist = 1e7, minLDcorr = 0.5,
                         overdispersion = 7)
hist(RD_09F2$alleleFreq, breaks = 50)
abline(v = 0.25, col = "blue")
abline(v = 0.5, col = "blue")

mat_09F2 <- GetWeightedMeanGenotypes(RD_09F2, omit1allelePerLocus = FALSE)
sum(is.na(mat_09F2))

# filter out the bulk of alleles that are not really segregating in the pop.
mat_09F2 <- mat_09F2[, RD_09F2$alleleFreq > 0.05 & RD_09F2$alleleFreq < 0.95]
dim(mat_09F2) # 5382 alleles
hist(colMeans(mat_09F2), breaks = 50) # many freq around 0.25, 0.5, not 0.75

# save(mat_09F2, file = file.path(outdir, "mat_09F2_HWE.RData"))

rm(RD_09F2)

load(file.path(outdir, "mat_09F2_HWE.RData"))

# Import Msa to split into 2x and 4x ####

RD_Msa <- readTagDigger(grep("Msa_counts", counts_files, value = TRUE),
                         dbfile = tagdb_file, possiblePloidies = list(2),
                         dbChrCol = "Chromosome", dbPosCol = "Position")
RD_Msa <- AddPCA(RD_Msa)
plot(RD_Msa)

accessions <- read.csv("C:/Users/lvclark/Documents/DOE study Msa/Seq/GBSv2_180110/all_accession_names.csv",
                       row.names = 1, stringsAsFactors = FALSE)
ploidies <- accessions[GetTaxa(RD_Msa), "Ploidy"]
names(ploidies) <- GetTaxa(RD_Msa)

table(ploidies)
colkey = c("black", "black", "red", "green", "blue")
names(colkey) <- c("", "unknown", "2x", "3x", "4x")

plot(RD_Msa, col = colkey[ploidies])
RD_Msa$PCA[ploidies %in% c("", "unknown"), 1:2]
GetTaxa(RD_Msa)[which(ploidies == "4x" & RD_Msa$PCA[,1] > 5)] # Golf Course 4x

ploidies[c("JY103", "JY103-HU", "JY263", "JY263-HU", "JY295", "JY298")] <- "2x"
ploidies[c("JY150")] <- "4x"

# get diploids, eliminating Msi
RD_Msa_2x <- SubsetByTaxon(RD_Msa, GetTaxa(RD_Msa)[ploidies == "2x" & RD_Msa$PCA[,2] < 20])

# tetraploids
RD_Msa_4x <- SubsetByTaxon(RD_Msa, GetTaxa(RD_Msa)[ploidies == "4x"])
RD_Msa_4x$possiblePloidies <- list(4L)

RD_Msa_2x # 360 taxa
RD_Msa_4x # 269 taxa

rm(RD_Msa)

# Msa diploids ####

RD_Msa_2x_prelim <- IterateHWE(RD_Msa_2x)

# odtest <- TestOverdispersion(RD_Msa_2x_prelim, to_test = 6:10)
# tiff(file.path(outdir, "Msa2xoverdispersion6.tiff"), compression = "lzw")
# qq(odtest[["6"]])
# dev.off()
# tiff(file.path(outdir, "Msa2xoverdispersion7.tiff"), compression = "lzw")
# qq(odtest[["7"]])
# dev.off()
# tiff(file.path(outdir, "Msa2xoverdispersion8.tiff"), compression = "lzw")
# qq(odtest[["8"]])
# dev.off()
# tiff(file.path(outdir, "Msa2xoverdispersion9.tiff"), compression = "lzw")
# qq(odtest[["9"]])
# dev.off()
# tiff(file.path(outdir, "Msa2xoverdispersion10.tiff"), compression = "lzw")
# qq(odtest[["10"]])
# dev.off()

rm(RD_Msa_2x_prelim)

RD_Msa_2x <- IteratePopStruct(RD_Msa_2x, overdispersion = 8)

mat_Msa_2x <- GetWeightedMeanGenotypes(RD_Msa_2x, omit1allelePerLocus = FALSE)

# get the same alleles that were retained from 09F2
mean(colnames(mat_09F2) %in% colnames(mat_Msa_2x)) # 90%
commnal <- colnames(mat_09F2)[colnames(mat_09F2) %in% colnames(mat_Msa_2x)]
mat_Msa2x <- matrix(0, nrow = nrow(mat_Msa_2x), ncol = ncol(mat_09F2),
                    dimnames = list(GetTaxa(RD_Msa_2x), colnames(mat_09F2)))
mat_Msa2x[,commnal] <- mat_Msa_2x[,commnal]

hist(colMeans(mat_Msa2x), breaks = 50) # broader range of al freq for smaller set
hist(colMeans(mat_Msa_2x), breaks = 50)

#save(mat_Msa2x, file = file.path(outdir, "mat_Msa2x.RData"))

rm(RD_Msa_2x)

load(file.path(outdir, "mat_Msa2x.RData"))

# Msa tetraploids ####

RD_Msa_4x_prelim <- IterateHWE(RD_Msa_4x)

odtest <- TestOverdispersion(RD_Msa_4x_prelim, to_test = 6:10)

# tiff(file.path(outdir, "Msa4xoverdispersion6.tiff"), compression = "lzw")
# qq(odtest[["6"]])
# dev.off()
# tiff(file.path(outdir, "Msa4xoverdispersion7.tiff"), compression = "lzw")
# qq(odtest[["7"]])
# dev.off()
# tiff(file.path(outdir, "Msa4xoverdispersion8.tiff"), compression = "lzw")
# qq(odtest[["8"]])
# dev.off()
# tiff(file.path(outdir, "Msa4xoverdispersion9.tiff"), compression = "lzw")
# qq(odtest[["9"]])
# dev.off()
# tiff(file.path(outdir, "Msa4xoverdispersion10.tiff"), compression = "lzw")
# qq(odtest[["10"]])
# dev.off()

rm(RD_Msa_4x_prelim)

# cull loci to save a little memory
keeploc <- unique(gsub("_[ACGT]*$", "", colnames(mat_09F2))) # 2008 loci
RD_Msa_4x <- SubsetByLocus(RD_Msa_4x, keeploc)

# run pipeline
RD_Msa_4x <- IteratePopStruct(RD_Msa_4x, overdispersion = 7)

mat_Msa_4x <- GetWeightedMeanGenotypes(RD_Msa_4x, omit1allelePerLocus = FALSE)

mat_Msa4x <- matrix(0, nrow = nrow(mat_Msa_4x), ncol = ncol(mat_09F2),
                    dimnames = list(GetTaxa(RD_Msa_4x), colnames(mat_09F2)))
mat_Msa4x[,commnal] <- mat_Msa_4x[,commnal]

hist(colMeans(mat_Msa4x), breaks = 50)

#save(mat_Msa4x, file = file.path(outdir, "mat_Msa4x.RData"))

load(file = file.path(outdir, "mat_Msa4x.RData"))

# Msi population ####

RD_Msi <- readTagDigger(grep("Msi_counts", counts_files, value = TRUE),
                        dbfile = tagdb_file, possiblePloidies = list(2),
                        dbChrCol = "Chromosome", dbPosCol = "Position")
RD_Msi <- SubsetByLocus(RD_Msi, keeploc) # cull unused loci
# eliminate blank and DH lines
RD_Msi <- 
  SubsetByTaxon(RD_Msi, GetTaxa(RD_Msi)[!GetTaxa(RD_Msi) %in% c("blank", "p196-150A-c", 
                                                                "p877-348-b", "IGR-2011-001")])
# merge taxa that were separate due to typo
RD_Msi <- MergeTaxaDepth(RD_Msi, c("PMS-504", "PMS-504-Msa"))

# test overdispersion
RD_Msi_prelim <- IterateHWE(RD_Msi)
odtest <- TestOverdispersion(RD_Msi_prelim, to_test = 6:10)
odtest2 <- TestOverdispersion(RD_Msi_prelim, to_test = 11:13)

# tiff(file.path(outdir, "Msi_overdispersion6.tiff"), compression = "lzw")
# qq(odtest[["6"]])
# dev.off()
# tiff(file.path(outdir, "Msi_overdispersion7.tiff"), compression = "lzw")
# qq(odtest[["7"]])
# dev.off()
# tiff(file.path(outdir, "Msi_overdispersion8.tiff"), compression = "lzw")
# qq(odtest[["8"]])
# dev.off()
# tiff(file.path(outdir, "Msi_overdispersion9.tiff"), compression = "lzw")
# qq(odtest[["9"]])
# dev.off()
# tiff(file.path(outdir, "Msi_overdispersion10.tiff"), compression = "lzw")
# qq(odtest[["10"]])
# dev.off()
# tiff(file.path(outdir, "Msi_overdispersion11.tiff"), compression = "lzw")
# qq(odtest2[["11"]])
# dev.off()
# tiff(file.path(outdir, "Msi_overdispersion12.tiff"), compression = "lzw")
# qq(odtest2[["12"]])
# dev.off()
# tiff(file.path(outdir, "Msi_overdispersion13.tiff"), compression = "lzw")
# qq(odtest2[["13"]])
# dev.off()

rm(RD_Msi_prelim)

RD_Msi <- IteratePopStruct(RD_Msi, overdispersion = 10)

mat_Msi_big <- GetWeightedMeanGenotypes(RD_Msi, omit1allelePerLocus = FALSE)

# get markers common to 09F2
commnal <- colnames(mat_09F2)[colnames(mat_09F2) %in% colnames(mat_Msi_big)]

mat_Msi <- matrix(0, nrow = nrow(mat_Msi_big), ncol = ncol(mat_09F2),
                    dimnames = list(GetTaxa(RD_Msi), colnames(mat_09F2)))
mat_Msi[,commnal] <- mat_Msi_big[,commnal]

hist(colMeans(mat_Msi), breaks = 50)

rm(RD_Msi)

#save(mat_Msi, file = file.path(outdir, "mat_Msi.RData"))

# Combined dataset ####
# uses the version of 09F2 that came from HWE pipeline
mat_combined <- rbind(mat_09F2, mat_Msa2x, mat_Msa4x, mat_Msi)

# visualize and explore dataset
mypca <- prcomp(mat_combined, scale. = TRUE)
mycol <- c(rep("purple", nrow(mat_09F2)), rep("blue", nrow(mat_Msa2x)),
           rep("green", nrow(mat_Msa4x)), rep("red", nrow(mat_Msi)))
plot(mypca$x[,1], mypca$x[,2], col = mycol)

hist(mat_combined)
hist(colMeans(mat_combined))
hist(colMeans(mat_combined[-(1:nrow(mat_09F2)),]))
sum(colMeans(mat_combined[-(1:nrow(mat_09F2)),]) == 0) # 242 alleles in 09F2 but not diversity panels

# alleles to toss for not being in diversity panels
notdiverse <- which(colMeans(mat_combined[-(1:nrow(mat_09F2)),]) == 0 |
                      colMeans(mat_combined[-(1:nrow(mat_09F2)),]) == 1)
str(notdiverse) #242
names(notdiverse)

mat_combined <- mat_combined[,-notdiverse]
dim(mat_combined) # 5140 markers remain

# write.csv(mat_combined, file = file.path(outdir, "combined_genotypes.csv"))
# write.csv(data.frame(Accession = rownames(mat_combined),
#                      Pop = c(rep("09F2", nrow(mat_09F2)), rep("Msa_2x", nrow(mat_Msa2x)),
#                                rep("Msa_4x", nrow(mat_Msa4x)), rep("Msi", nrow(mat_Msi)))),
#           file = file.path(outdir, "accessions_populations.csv"),
#           row.names = FALSE)

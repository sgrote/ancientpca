#' VCF file to genotype matrix
#' Load vcf file into R and extract a genotype matrix coded with 0,1,2
#' @param vcf_file 
#' @param max_missing_snp Set a max. percentage of missingness per SNP allowed. Default is set to 1 (i.e. 100%)
#' @param max_missing_sample Set a max. percentage of missingness per sample allowed.  Default is set to 1 (i.e. 100%)
#' @export

vcfToGenotypeMatrix <- function(vcf_file, max_missing_snp=1, max_missing_sample=1){
	vcf <- read.vcfR( vcf_file )
	head(vcf)

	# extract genotype matrix 
	gt <- extract.gt(vcf)

	#recode 0,1,2
	gt[gt == "0/0"] <- 0
	gt[gt == "0/1"] <- 1
	gt[gt == "1/0"] <- 1
	gt[gt == "1/1"] <- 2
	gt <- t(gt)
	storage.mode(gt) <- "numeric"

	# remove samples (rows), that have a greater that the set max. missingness
	if (max_missing_sample < 1){
		print("removing individuals:")
		print(rownames(as.data.frame(which(rowMeans(is.na(gt))  > max_missing_sample))))
		gt <- gt[-which(rowMeans(is.na(gt)) >= max_missing_sample), ]
	}
	# remove snps (columns), that have a greater that the set max. missingness
	if (max_missing_snp < 1){
		gt <- gt[, -which(colMeans(is.na(gt)) >= max_missing_snp)]
	}

	str(gt)

	return(gt)
}


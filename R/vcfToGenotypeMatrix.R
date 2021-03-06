#' VCF file to genotype matrix
#'
#' Load vcf file into R and extract a genotype matrix (sample X SNP) coded with 0,1,2,<NA>. 
#' @param vcf_file VCF file with genotypes only (i.e. 0/0, 0/1, etc.). Missing values can be coded as "." or "./."
#' @param max_missing_sample Set a max. percentage of missingness per sample allowed.  Default is set to 1 (i.e. 100 percent). Samples with more missing genotypes will be removed.
#' @param max_missing_snp Set a max. percentage of missingness per SNP allowed. Default is set to 1 (i.e. 100 percent). Sites where more than `max_missing_snp` of the samples have no genotypes, will be removed. This is done after samples with too many missing data are removed.
#' @return matrix (sample x snp) coded 0,1,2, or NA
#' @import utils
#' @export

vcfToGenotypeMatrix <- function(vcf_file, max_missing_snp=1, max_missing_sample=1){

	vcf <- vcfR::read.vcfR( vcf_file )
	head(vcf)

	# extract genotype matrix
	gt <- vcfR::extract.gt(vcf)

	#recode 0,1,2
	gt[gt == "0/0"] <- 0
	gt[gt == "0/1"] <- 1
	gt[gt == "1/0"] <- 1
	gt[gt == "1/1"] <- 2
	gt <- t(gt)
	storage.mode(gt) <- "numeric"

	# remove samples (rows), that have a greater than the set max. missingness
	if (max_missing_sample < 1){
		print("Removing individuals:")
		print(rownames(as.data.frame(which(rowMeans(is.na(gt))  >= max_missing_sample))))
		gt <- gt[-which(rowMeans(is.na(gt)) > max_missing_sample), ]
	}

	# remove snps (columns), that have a greater than the set max. missingness
	if (max_missing_snp < 1){
		gt <- gt[, -which(colMeans(is.na(gt)) > max_missing_snp)]
	}

	print("Genotype Matrix dimensions:")
	print(paste0(dim(gt)[1], " samples x ", dim(gt)[2], " SNPs"))
	
	# save thresholds
	attr(gt, "max_missing_snp") <- max_missing_snp
	attr(gt, "max_missing_sample") <- max_missing_sample
	
	return(gt)
}

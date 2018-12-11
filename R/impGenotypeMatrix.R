#' Imputed missing genotypes in incomplete genotype matrix
#'
#' Returns a matrix with imputed missing data in a 0,1,2 coded genotype matrix
#' @param genotype_matrix Genotype matrix (sample x snps) coded with 0,1,2, or NA
#' @export

impGenotypeMatrix <- function(genotype_matrix){

	library(softImpute)
	fits <- softImpute(genotype_matrix, trace=TRUE, type="svd")
	gt_imputed <- complete(genotype_matrix, fits)


	return(gt_imputed)
}

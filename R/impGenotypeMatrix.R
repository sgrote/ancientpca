#' Impute missing genotypes in incomplete genotype matrix
#'
#' Returns a matrix with imputed missing data in a 0,1,2 coded genotype matrix. Non-variance columns are automatically removed.
#' @param genotype_matrix Genotype matrix (sample x snps) coded with 0,1,2, or <NA>
#' @return matrix (sample x snp) with genotype "probabilities" ranging from 0 to 2
#' @importFrom stats var
#' @export

impGenotypeMatrix <- function(genotype_matrix){

	fits <- softImpute::softImpute(genotype_matrix, trace=TRUE, type="svd")
	gt_imputed <- softImpute::complete(genotype_matrix, fits)


	# remove non-variance columns - only needed if snps or samples were removed
	gt_imputed <- gt_imputed[ , (apply(gt_imputed, 2, var) != 0)]


	return(gt_imputed)
}

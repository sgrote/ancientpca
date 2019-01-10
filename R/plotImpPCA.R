#' Plot PCAs 1-6
#'
#' Exports plots of PCs 1~2, 3~4, 5~6 to PDF file
#' @param pca_obj PCA object from R funtion prcomp()
#' @param original_matrix Original genotype matrix coded 0,1,2, or <NA>
#' @param meta_file CSV file with two columns: sample-ID and population. Might contain additional columns.
#' @param output_pca_pdf Output file with ".pdf"
#' @import ggplot2
#' @import utils
#' @export


plotImpPCA <- function(pca_obj, original_matrix, meta_file, output_pca_pdf){
    
	# extract rotated loadings
	loads <- pca_obj$x

	# add meta data
	meta <- read.table(meta_file, sep = ",", header=TRUE)
	miss_meta <- rownames(loads)[!(rownames(loads) %in% meta[,1])]
	if (length(miss_meta) > 0){
	    stop("missing metadata for samples: ", paste(miss_meta, collapse=", "))
	}
	meta_subset <- meta[match(rownames(loads), meta[,1]), ]

	#get missing stats on samples in original data
	density <- character(0)
	for (i in 1:length(meta_subset[,1])) density[i] <- round(1 -(sum(is.na(original_matrix[i,]))/dim(original_matrix)[2]), 1)

	# make a data.frame with individuals and PC1 and PC2
	tab <- data.frame(sample.id = meta_subset[,1],
    	population = meta_subset[,2],
    	density = density,
    	EV1 = loads[,1],    # the first eigenvector
    	EV2 = loads[,2],
    	EV3 = loads[,3],
    	EV4 = loads[,4],
    	EV5 = loads[,5],
    	EV6 = loads[,6],
    	stringsAsFactors = FALSE)
	head(tab)


	pc.percent <- pca_obj$sdev^2/sum(pca_obj$sdev^2)*100


	# Draw
	# PC 1 and 2

	plot1.0 <- ggplot(tab, aes(EV1, EV2,
                         color= tab$population,
                         label = tab$sample.id, size=tab$density, group=tab$sample.id)) +
    	ggrepel::geom_label_repel(aes(EV1, EV2, label = tab$sample.id)) +
    	scale_shape_manual(values=1:8) +
    	geom_point(size=1, stroke = 2) +
    	xlab(paste("PC 1 (", round(pc.percent, 2)[1], "%)", sep = "")) +
    	ylab(paste("PC 2 (", round(pc.percent, 2)[2], "%)", sep = "")) +
    	theme_bw() +
    	theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(colour = "black", size=0.5)) +
    	theme(legend.box = "horizontal") +
    	labs(col="population/site", size="density of data") +
    	theme(legend.text=element_text(size=13))

	plot1 <- plot1.0 + theme(legend.position = "none")


	# PC 3 and 4
	plot2 <- ggplot(tab, aes(EV3, EV4,
                         color= tab$population,
                         label = tab$sample.id, size=tab$density)) +
    	ggrepel::geom_label_repel(aes(EV3, EV4, label = tab$sample.id)) +
    	scale_shape_manual(values=1:8) +
    	geom_point(size=1, stroke = 2) +
    	xlab(paste("PC 3 (", round(pc.percent, 2)[3], "%)", sep = "")) +
    	ylab(paste("PC 4 (", round(pc.percent, 2)[4], "%)", sep = "")) +
    	theme_bw() +
    	theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(colour = "black", size=0.5)) +
    	theme(legend.position = "none")


	# PC 4 and 5
	plot3 <- ggplot(tab, aes(EV5, EV6,
                         color= tab$population,
                         label = tab$sample.id, size=tab$density)) +
    	ggrepel::geom_label_repel(aes(EV5, EV6, label = tab$sample.id)) +
    	scale_shape_manual(values=1:8) +
    	geom_point(size=1, stroke = 2) +
    	xlab(paste("PC 5 (", round(pc.percent, 2)[5], "%)", sep = "")) +
    	ylab(paste("PC 6 (", round(pc.percent, 2)[6], "%)", sep = "")) +
    	theme_bw() +
    	theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(colour = "black", size=0.5)) +
    	theme(legend.position = "none")

	# legend
	legend <- cowplot::get_legend(plot1.0)

	# combine plots
	p <- cowplot::plot_grid(plot1, plot2, plot3, legend, labels = c("A", "B", "C"))

	# now add the title
	title_string = paste0("max. ", attr(original_matrix, "max_missing_snp")*100, "% missingness per SNP and max. ", attr(original_matrix, "max_missing_sample")*100, "% missingness per sample")
	title <- cowplot::ggdraw() + cowplot::draw_label(title_string, fontface = 'bold')

	p2 <- cowplot::plot_grid(title, p, ncol = 1, rel_heights = c(0.04, 1)) # rel_heights values control title margins

	# save plots to one PDF using "cowplot"
	cowplot::save_plot(output_pca_pdf, p2, ncol = 3, nrow=3)
}

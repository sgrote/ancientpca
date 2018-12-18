#' Plot PCAs 1-6
#'
#' Exports plots of PCs 1~2, 3~4, 5~6 to PDF file
#' @param max_missing_snp Set a max. percentage of missingness per SNP allowed. Default is set to 1 (i.e. 100 percent)
#' @param max_missing_sample Set a max. percentage of missingness per sample allowed.  Default is set to 1 (i.e. 100 percent)
#' @param meta_file CVS file that contains headers "Sample_ID" and "Population". Can contain other columns as well
#' @export


plotImpPCA <- function(pca_obj, imputed_matrix, original_matrix, max_missing_snp=1, max_missing_sample=1, meta_file, output_pca_filename){
	library(ggplot2)
	library(ggrepel)
	library(cowplot)

	# add meta data

	meta <- read.table(meta_file, sep = ",", header=TRUE)
	meta_subset <- meta[meta$Sample_ID %in% rownames(gt_imputed), ]


	#get missing stats on samples in original data
	density <- character(0)
	for (i in 1:length(meta_subset$Sample_ID)) density[i] <- round(1 -(sum(is.na(gt[i,]))/dim(gt)[2]), 1)

	# make a data.frame with individuals and PC1 and PC2
	tab <- data.frame(sample.id = meta_subset$Sample_ID,
    	population = meta_subset$Population,
    	density = density,
    	EV1 = pca_gt$x[,1],    # the first eigenvector
    	EV2 = pca_gt$x[,2],
    	EV3 = pca_gt$x[,3],
    	EV4 = pca_gt$x[,4],
    	EV5 = pca_gt$x[,5],
    	EV6 = pca_gt$x[,6],
    	stringsAsFactors = FALSE)
	head(tab)


	pc.percent <- pca_gt$sdev^2/sum(pca_gt$sdev^2)*100


	# Draw
	# PC 1 and 2

	plot1.0 <- ggplot(tab, aes(EV1, EV2,
                         color= tab$population,
                         label = tab$sample.id, size=tab$density, group=tab$sample.id)) +
    	geom_label_repel(aes(EV1, EV2, label = tab$sample.id)) +
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
    	geom_label_repel(aes(EV3, EV4, label = tab$sample.id)) +
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
    	geom_label_repel(aes(EV5, EV6, label = tab$sample.id)) +
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
	legend <- get_legend(plot1.0)

	# combine plots
	p <- plot_grid(plot1, plot2, plot3, legend, labels = c("A", "B", "C"))

	# now add the title
	title_string = paste0("max. ", toString(max_missing_snp*100), "% missingness per SNP and max. ", toString(max_missing_sample*100), "% missingness per sample")
	title <- ggdraw() + draw_label(title_string, fontface = 'bold')

	p2 <- plot_grid(title, p, ncol = 1, rel_heights = c(0.04, 1)) # rel_heights values control title margins

	# save plots to one PDF using "cowplot"
	save_plot(output_pca_filename, p2, ncol = 3, nrow=3)
}

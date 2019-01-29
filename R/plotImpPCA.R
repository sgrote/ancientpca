

#' Plot 1 PCA plot
#' @param loads matrix of PC loads
#' @param pca_x PC nr for x-axis
#' @param pca_y PC nr for y-axis
#' @param meta metadata for size and legend
#' @param pc_percent proportion of variance explained for all PCs
#' @return PCA plot

pca_plot <- function(loads, pca_x, pca_y, meta, pc_percent){
    plotty <- ggplot(mapping = aes(loads[,pca_x], loads[,pca_y], color=meta$population,
	label=meta$sample_id, size=meta$dens, group=meta$sample_id)) +
	ggrepel::geom_label_repel(aes(loads[,pca_x], loads[,pca_y], label=meta$sample_id)) +
	scale_size_continuous(range=c(3,6)) +
	geom_point(size=1, stroke=2) +
	xlab(paste("PC ", pca_x, " (" , round(pc_percent[pca_x],2), "%)", sep = "")) +
	ylab(paste("PC ", pca_y, " (", round(pc_percent[pca_y],2), "%)", sep = "")) +
	theme_bw() +
	theme(panel.grid.major = element_blank(),
	    panel.grid.minor = element_blank(),
	    panel.background = element_rect(colour = "black", size=0.5)) +
	theme(legend.box = "horizontal") +
	labs(col="population/site", size="density of data (SNPs covered)") +
	theme(legend.text=element_text(size=13))
    return(plotty)
}



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
    pc_percent <- pca_obj$sdev^2/sum(pca_obj$sdev^2)*100

    # add meta data
    meta <- read.table(meta_file, sep = ",", header=TRUE)
    miss_meta <- rownames(loads)[!(rownames(loads) %in% meta[,1])]
    if (length(miss_meta) > 0){
	stop("missing metadata for samples: ", paste(miss_meta, collapse=", "))
    }
    meta <- meta[match(rownames(loads), meta[,1]), ]

    # get missing stats on samples in original data (already filtered with max_miss*)
    dens <- 1 - (rowMeans(is.na(original_matrix)))

    # make a data.frame with individuals and PCs
    plot_meta <- data.frame(meta[,1], meta[,2], dens, stringsAsFactors=FALSE)
    colnames(plot_meta) <- c("sample_id", "population", "dens")
    
    # PCA plots for different PCs
    plot1 <- pca_plot(loads, 1, 2, plot_meta, pc_percent)
    plot2 <- pca_plot(loads, 3, 4, plot_meta, pc_percent) + theme(legend.position = "none")
    plot3 <- pca_plot(loads, 5, 6, plot_meta, pc_percent) + theme(legend.position = "none")

    # legend
    legend <- cowplot::get_legend(plot1)
    plot1 <- plot1 + theme(legend.position = "none")

    # combine plots
    p <- cowplot::plot_grid(plot1, plot2, plot3, legend, labels = c("A", "B", "C"))

    # now add the title
    title_string = paste0("max. ", attr(original_matrix, "max_missing_snp")*100, "% missingness per SNP and max. ", attr(original_matrix, "max_missing_sample")*100, "% missingness per sample")
    title <- cowplot::ggdraw() + cowplot::draw_label(title_string, fontface = 'bold')

    # rel_heights values control title margins
    p2 <- cowplot::plot_grid(title, p, ncol=1, rel_heights=c(0.04,1)) 

    # save plots to one PDF using "cowplot"
    cowplot::save_plot(output_pca_pdf, p2, ncol=3, nrow=3)
    
    # cowplot in addition opens an empty plotting device
    dev.off()
    
    return(invisible(NULL))
}

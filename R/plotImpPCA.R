

#' Plot 1 PCA plot
#' @param loads matrix of PC loads
#' @param pc_x PC nr for x-axis
#' @param pc_y PC nr for y-axis
#' @param meta metadata for size and legend
#' @param pc_percent proportion of variance explained for all PCs
#' @return PCA plot

pca_plot <- function(loads, pc_x, pc_y, meta, pc_percent){
    plotty <- ggplot(mapping = aes(loads[,pc_x], loads[,pc_y], color=meta$population,
	label=meta$sample_id, size=meta$dens, group=meta$sample_id)) +
	ggrepel::geom_label_repel(aes(loads[,pc_x], loads[,pc_y], label=meta$sample_id)) +
	scale_size_continuous(range=c(2,4.5)) +
	geom_point(size=1, stroke=2) +
	xlab(paste("PC ", pc_x, " (", round(pc_percent[pc_x],2), "%)", sep = "")) +
	ylab(paste("PC ", pc_y, " (", round(pc_percent[pc_y],2), "%)", sep = "")) +
	theme_bw() +
	theme(aspect.ratio=1) +
	coord_fixed() +
	theme(panel.grid.major = element_blank(),
	    panel.grid.minor = element_blank(),
	    panel.background = element_rect(colour = "black", size=0.5)) +
	theme(legend.box = "horizontal") +
	labs(col="population/site", size="density of data") +
	theme(legend.text=element_text(size=13))
    return(plotty)
}



#' Plot PCAs
#' Exports combined plot to PDF file
#' @param pca_obj PCA object from R funtion prcomp()
#' @param original_matrix Original genotype matrix coded 0,1,2, or <NA>
#' @param meta_file CSV file with two columns: sample-ID and population. Might contain additional columns.
#' @param output_pca_pdf Output file with ".pdf"
#' @param pc_max Even integer, last PC to plot, e.g. if pc_max=4 it will plot PC1 vs PC2 and PC3 vs PC4.
#' @import ggplot2
#' @import utils
#' @export

plotImpPCA <- function(pca_obj, original_matrix, meta_file, output_pca_pdf, pc_max=6){
    
    # extract rotated loadings
    loads <- pca_obj$x
    
    # check that pc_max is an even number and <= nr. of PCs 
    val_pc_max = seq(2, ncol(loads), 2)
    if (! pc_max %in% val_pc_max){
	stop("'pc_max' has to be one of {", paste(val_pc_max, collapse=", "), "}.")
    }
    
    # explained variance per PC
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
    plot_list <- list()
    for (i in 1:(pc_max/2)){
	pc_x <- i * 2 - 1
	pc_y <- i * 2
	plot_list[[i]] <- pca_plot(loads, pc_x, pc_y, plot_meta, pc_percent)
	# keep legend in first plot to extract it below
	if (i > 1){
	    plot_list[[i]] <- plot_list[[i]] + theme(legend.position = "none")
	}
    }
    
    # extract and add legend in separate plot
    legend <- cowplot::get_legend(plot_list[[1]])
    plot_list[[1]] <- plot_list[[1]] + theme(legend.position = "none")
    plot_list <- c(plot_list, list(legend))

    # combine plots
    labels <- LETTERS[1:(pc_max/2)]
    n <- length(plot_list)
    ncol <- ceiling(sqrt(n))
    nrow <- ceiling(n/ncol)
    p1 <- cowplot::plot_grid(plotlist=plot_list, labels=labels, nrow=nrow, ncol=ncol)
    
    # now add the title
    title_string = paste0("max. ", attr(original_matrix, "max_missing_snp")*100, "% missingness per SNP and max. ", attr(original_matrix, "max_missing_sample")*100, "% missingness per sample")
    title <- cowplot::ggdraw() + cowplot::draw_label(title_string, fontface='bold')
    # rel_heights values control title margins
    p2 <- cowplot::plot_grid(title, p1, ncol=1, rel_heights=c(0.1/nrow,1))

    # save plots to one PDF using "cowplot"
    cowplot::save_plot(output_pca_pdf, p2, base_width=5, base_height=5, nrow=nrow, ncol=ncol)
    
    # cowplot in addition opens an empty plotting device
    dev.off()
    
    return(invisible(NULL))
}


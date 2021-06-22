
#---- GSEA functions

#Written by Monica Hannani - https://github.com/monicahannani/CovidEpiMap/blob/main/sc_source/sc_source.R

run_gsea = function(bg.genes, stats, category, plot.title = NULL, subcategory = NULL, out.dir = '.', file.prefix, n = 30){
  suppressPackageStartupMessages(library(msigdbr))
  suppressPackageStartupMessages(library(fgsea))
  suppressPackageStartupMessages(library(dplyr))
  
  # Fetch geneset
  geneSets = msigdbr(species = 'Homo sapiens', category = category, subcategory = subcategory)
  geneSets = geneSets[geneSets$human_gene_symbol %in% bg.genes,]
  m_list = geneSets %>% split(x = .$human_gene_symbol, f = .$gs_name)
  
  # Run GSEA
  gsea = fgsea(pathways = m_list, stats = stats, minSize = 10, eps = 0.0)
  order = order(gsea$padj, decreasing = FALSE)
  
  # Plot
  file.name = paste0(out.dir, '/', file.prefix, '_gsea_', paste0(c(category, subcategory), collapse = '_'))
  
  pdf(file = paste0(file.name, '.pdf'), width = 10, height = 9)
  print(plot_go(gsea.res = gsea,
                gsea.res.order = order, 
                n = n, 
                plot.title = plot.title))
  dev.off()
  
  # Write to file
  write.table(gsea[order, -8], 
              file = paste0(file.name, '.txt'), 
              sep = '\t', 
              row.names = FALSE, 
              quote = FALSE)
}

#---- Plot GO terms from GSEA

plot_go = function(gsea.res, gsea.res.order, plot.title = NULL, n = 20){
  suppressPackageStartupMessages(library(ggplot2))
  
  # Format GOs
  plot.table = head(gsea.res[gsea.res.order,], n = n)
  plot.table$pathway = sub('GO_', '', plot.table$pathway)
  plot.table$pathway = gsub('_', ' ', plot.table$pathway)
  
  p = ggplot(plot.table,
             aes(x = NES, y = pathway)) +
    geom_point(aes(size = NES, color = padj)) +
    theme_bw(base_size = 8) +
    ylab(NULL) +
    ggtitle(plot.title) +
    scale_colour_gradient2(low = 'red', 
                           mid = 'lightgrey', 
                           high = 'blue', 
                           midpoint = 0.05, 
                           limits = c(0,0.1), 
                           oob = scales::squish)
  
  return(p)
}

#---- Create expression scatter plot of gene lists

expr.plot <- function(cluster, genes, pval=TRUE) {
  c <- metadata %>% filter(seurat_clusters == cluster)
  
  c.df <- data@assays$RNA@data %>% as.data.frame()
  
  c.df <- c.df[,colnames(c.df) %in% row.names(c),]
  
  #Select genes of interest
  df.goi <- c.df[row.names(c.df) %in% genes,] %>% t()
  
  #merge metadata info back onto df
  c.merge <- merge(df.goi, c, all=TRUE, by="row.names")
  
  #Create short vector of only gene names that were found
  genes.short <- colnames(df.goi)
  
  #Put in long format
  c.merge <- c.merge %>% gather(Gene, Value, genes.short)
  
  if(pval){
    ggplot(c.merge, aes(VTbin, Value, colour=VTbin)) +
      geom_jitter() +
      facet_wrap(~Gene) +
      stat_compare_means(method = "t.test", label.y = 4) +
      ylim(0,5)
  } else {
    #Plot expression differences
    ggplot(c.merge, aes(VTbin, Value, colour=VTbin)) +
      geom_jitter() +
      facet_wrap(~Gene)
  }
}

#---- Create expression scatter plot of gene lists

score.comparison <- function(cluster){
  #Add senescence score
  score <- AddModuleScore(
    object = cluster,
    features = list(senescent.g),
    ctrl = 5,
    name = 'senes')
  
  #Subset metadata with senescence scores
  score.meta <- score@meta.data
  
  #Set pdf name
  pdfn <- paste("senescence/cluster", unique(score.meta$cluster), "_", "senescence_score.pdf", sep="")
  
  #Write plot
  pdf(file = pdfn, width = 10, height = 9)
  #Plot boxplot of senes results colouring by infection with BK
  print(ggplot(score.meta, aes(x=VTbin, y=senes1, fill=VTbin)) +
          geom_boxplot() +
          stat_compare_means(method = "t.test", label.y = 0.6) +
          ggtitle(paste("Cluster",  unique(score.meta$cluster),": Comparison of senescence scores between virally infected cells and controls", sep="")))
  dev.off()
  
}

#---- Create plot averaging expression over gene set
gsea.expr.plot <- function(cluster, genes, group) {
  c <- metadata %>% filter(seurat_clusters == cluster)
  
  c.df <- data@assays$RNA@data %>% as.data.frame()
  
  c.df <- c.df[,colnames(c.df) %in% row.names(c),]
  
  #Select genes of interest
  df.goi <- c.df[row.names(c.df) %in% genes,] %>% t()
  
  means <- rowMeans(df.goi, na.rm=TRUE) %>% as.data.frame()
  
  colnames(means) <- "avg"
  
  df.goi <- cbind(df.goi, means)
  
  #merge metadata info back onto df
  c.merge <- merge(df.goi, c, all=TRUE, by="row.names")
  
  #Plot expression differences
  ggplot(c.merge, aes_string(group, "avg", colour=group)) +
    geom_jitter() +
    stat_summary(fun = median, fun.min = median, fun.max = median,
                 geom = "crossbar", width = 0.5, colour="#616A6B")
}
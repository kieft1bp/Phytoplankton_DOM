library(data.table)
library(reshape2)
library(gridExtra)
library(dplyr)
library(ggplot2)
library(grid)
library(RColorBrewer)
library(vegan)
library(ape)
library(splitstackshape)
library(cowplot)
library(Ternary)
library(genefilter)
library(gdata)

# Main text figures
plot_fig1a = "T"     ## Per-sample enrv and lf summary
plot_fig1b = "T"     ## Unlabeled taxonomic composition PCoA plot
plot_fig1c = "T"     ## Labeled taxonomic composition PCoA plot
plot_fig2a = "T"     ## Labeled taxonomic composition barplot
plot_fig2b = "T"     ## Per-lineage per-treatment ENRV histograms
plot_fig3a = "T"     ## Proteome function ternary plots of controls
plot_fig3b = "T"     ## Proteome function ternary plots of treatments
plot_fig4a = "T"     ## Assimilation summary across treatments (ENRV*LF)
plot_fig4b = "T"     ## DEMIC Growth rate estimations from binned populations 
plot_fig4c = "T"     ## Estimated C:N requirements from genes and proteins per sig enriched genus

# Supplemental figures
plot_sfig1 = "T"     ## Unlabeled taxonomic composition PCoA plot with controls
plot_sfig2 = "T"     ## Significantly-enriched genera Z-scores of ENRV
plot_sfig4 = "T"     ## Unlabeled taxonomic composition barplot
plot_sfig6 = "T"     ## Functional deviation scatterplots of taxonomic Orders
plot_sfig7 = "T"     ## Per-lineage boxplots of proteome functional composition

#####################################
############# FIGURE 1a ############# 
#####################################
if(plot_fig1a == "T"){
  
  # Read in data
  all_data = read.csv("all_protein_data.csv", header=T, sep=",", stringsAsFactors = F, comment.char = "")
  
  # Remove extraneous information
  sample_by_enrichment_df = aggregate(. ~ sample, data=all_data[,c(1,12:ncol(all_data))], FUN=sum)
  sample_by_enrichment_df = sample_by_enrichment_df[!grepl("^t",sample_by_enrichment_df$sample),]
  
  # Define labeled (>= 2%) and unlabeled (< 2%) fractions 
  samples = sample_by_enrichment_df[,1]
  all_proteins = sample_by_enrichment_df[,2:ncol(sample_by_enrichment_df)]
  unlabeled_proteins = sample_by_enrichment_df[,2:3]
  labeled_proteins = sample_by_enrichment_df[,4:ncol(sample_by_enrichment_df)]
  
  # Calculate ENRV
  enrv_multiplier = as.numeric(2:100)
  enrv = data.frame(mapply(`*`,labeled_proteins,enrv_multiplier))
  enrv = rowSums(enrv)/rowSums(labeled_proteins)
  
  # Calculate LF
  lf = (rowSums(labeled_proteins)/rowSums(all_proteins))*100
  
  # Prepare data for plotting
  fig1a_data = data.frame(samples,enrv,lf)
  colnames(fig1a_data) = c("sample","enrv","lf")
  fig1a_data$sample = substr(fig1a_data$sample,1,nchar(as.vector(fig1a_data$sample))-1)
  
  # Calculate stats
  enrv_model <- lm(enrv ~ sample, data=fig1a_data)
  anova(enrv_model)
  enrv_anova = aov(enrv_model)
  TukeyHSD(enrv_anova)
  
  lf_model <- lm(lf ~ sample, data=fig1a_data)
  anova(lf_model)
  lf_anova = aov(lf_model)
  TukeyHSD(lf_anova)
  
  # Plot
  fig1a_data = melt(fig1a_data)
  fig1a_data$sample = factor(fig1a_data$sample,levels = c("DLy","DEx","CLy","CEx"))
  fig1a_plot = ggplot(data=fig1a_data) +
    geom_violin(aes(x=sample,y=value,linetype=variable,fill=sample),draw_quantiles = 0.5,scale="width",position = position_dodge(width = 0.8), show.legend = F) +
    #geom_boxplot(aes(x=sample,y=value,linetype=variable,fill=sample),position = position_dodge(width = 0.8), show.legend = F) +
    scale_linetype_manual(values = c(1,9)) +
    scale_fill_manual(values = c("gray20","gray40","gray60","gray80")) +
    theme_minimal() +
    scale_y_continuous(name="ENRV (top) and LF (bottom) %",limits = c(0,100),expand = c(0,0)) +
    theme(axis.title.x=element_blank(),
          axis.title.y=element_text(color="black",size=10),
          legend.title = element_blank(),
          strip.background = element_blank(),
          strip.text.x = element_blank(),
          plot.title = element_blank(),
          axis.line.y = element_line(color = "black"),
          axis.line.x = element_line(color = "black"),
          axis.ticks = element_line(color = "black"),
          axis.text.x = element_text(color = "black",size=12),
          axis.text.y = element_text(color = "black",size=12),
          axis.title = element_text(color = "black"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
  
  ggsave("fig1a.svg", fig1a_plot, device = "svg",scale = 1, width = 3, height = 2.5, units = "in", dpi = 600)
}

#####################################
############# FIGURE 1b ############# 
#####################################
if(plot_fig1b == "T"){
  
  # Read in data
  all_data = read.csv("all_protein_data.csv", header=T, sep=",", stringsAsFactors = F, comment.char = "")
  
  # Aggregate rows
  sample_by_enrichment_df = aggregate(. ~ species + sample, data=all_data[,c(1,11,12:ncol(all_data))], FUN=sum)
  sample_by_enrichment_df = sample_by_enrichment_df[!grepl("^t",sample_by_enrichment_df$sample),]
  
  # Aggregate columns
  sample_by_enrichment_df$sum = rowSums(sample_by_enrichment_df[,c(3:4)])
  unlabeled_proteins = sample_by_enrichment_df[,c(1:2,ncol(sample_by_enrichment_df))]

  # Unmelt data to create relative abundance matrix
  unlabeled_proteins_unmelted = dcast(data = unlabeled_proteins,formula = sample~species,fun.aggregate = sum,value.var = "sum")
  rownames(unlabeled_proteins_unmelted) = unlabeled_proteins_unmelted$sample
  unlabeled_proteins_unmelted$sample = NULL
  unlabeled_proteins_unmelted = as.data.frame(t(unlabeled_proteins_unmelted))
  cols = as.character(colnames(unlabeled_proteins_unmelted))
  rows = as.character(rownames(unlabeled_proteins_unmelted))
  unlabeled_proteins_unmelted = unlabeled_proteins_unmelted %>% mutate_all(as.character)
  unlabeled_proteins_unmelted = sapply(unlabeled_proteins_unmelted, as.numeric)
  unlabeled_proteins_unmelted = sweep(unlabeled_proteins_unmelted, 2, colSums(unlabeled_proteins_unmelted), '/')
  rownames(unlabeled_proteins_unmelted) = rows
  colnames(unlabeled_proteins_unmelted) = cols

  # Calculate distance matrix
  samples_dist = vegdist(t(unlabeled_proteins_unmelted), method = "bray")
  
  # Calculate PCoA on distance matrix
  pcoa = pcoa(samples_dist)
  axis_1_variance = round(pcoa$values$Relative_eig[1]*100,digits = 2)
  axis_2_variance = round(pcoa$values$Relative_eig[2]*100, digits = 2)
  pcoa_vectors = data.frame(pcoa$vectors[,1:2])
  vectors_to_plot = pcoa_vectors
  colnames(vectors_to_plot) = c("Axis.1","Axis.2")
  vectors_to_plot$sample = rownames(vectors_to_plot)
  vectors_to_plot$treatment = substr(vectors_to_plot$sample,1,nchar(as.vector(vectors_to_plot$sample))-1)

  fig1b_plot = ggplot(data=vectors_to_plot, aes(x=Axis.1, y = Axis.2,fill=treatment, color=treatment)) +
    geom_hline(yintercept = 0, lty="dotted")+
    geom_vline(xintercept = 0, lty="dotted")+
    geom_polygon(size=0.5,alpha=0.2) +
    geom_point(shape=21, size= 2.5) +
    scale_color_manual(name="",guide = guide_legend(reverse = TRUE),values = c("CEx"="gray80","CLy"="gray60","DEx"="gray40","DLy"="gray20")) +
    scale_fill_manual(name="",guide = guide_legend(reverse = TRUE),values = c("CEx"="gray80","CLy"="gray60","DEx"="gray40","DLy"="gray20")) +
    xlab(paste("PCoA 1 (",axis_1_variance,"%)",sep="")) + 
    ylab(paste("PCoA 2 (",axis_2_variance,"%)",sep="")) + 
    xlim(c(-0.3,0.4)) +
    ylim(c(-0.2,0.14)) +
    theme(panel.border = element_rect(colour = "black", fill=NA),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.title = element_text(size=12,color="black"),
          axis.text.x=element_text(size=12,color="black"), 
          axis.text.y=element_text(size=12,color="black"))
  
  ggsave("fig1b.svg", fig1b_plot, device = "svg",scale = 1, width = 3.5, height = 2.5, units = "in", dpi = 600)
  
}
  
#####################################
############# FIGURE 1c ############# 
#####################################
if(plot_fig1c == "T"){
  
  # Read in data
  all_data = read.csv("all_protein_data.csv", header=T, sep=",", stringsAsFactors = F, comment.char = "")
  
  # Aggregate rows
  sample_by_enrichment_df = aggregate(. ~ species + sample, data=all_data[,c(1,11,12:ncol(all_data))], FUN=sum)
  sample_by_enrichment_df = sample_by_enrichment_df[!grepl("^t",sample_by_enrichment_df$sample),]
  
  # Aggregate columns
  sample_by_enrichment_df$sum = rowSums(sample_by_enrichment_df[,c(5:ncol(sample_by_enrichment_df))])
  labeled_proteins = sample_by_enrichment_df[,c(1:2,ncol(sample_by_enrichment_df))]
  
  # Unmelt data to create relative abundance matrix
  unlabeled_proteins_unmelted = dcast(data = labeled_proteins,formula = sample~species,fun.aggregate = sum,value.var = "sum")
  rownames(unlabeled_proteins_unmelted) = unlabeled_proteins_unmelted$sample
  unlabeled_proteins_unmelted$sample = NULL
  unlabeled_proteins_unmelted = as.data.frame(t(unlabeled_proteins_unmelted))
  cols = as.character(colnames(unlabeled_proteins_unmelted))
  rows = as.character(rownames(unlabeled_proteins_unmelted))
  unlabeled_proteins_unmelted = unlabeled_proteins_unmelted %>% mutate_all(as.character)
  unlabeled_proteins_unmelted = sapply(unlabeled_proteins_unmelted, as.numeric)
  unlabeled_proteins_unmelted = sweep(unlabeled_proteins_unmelted, 2, colSums(unlabeled_proteins_unmelted), '/')
  rownames(unlabeled_proteins_unmelted) = rows
  colnames(unlabeled_proteins_unmelted) = cols
  
  # Calculate distance matrix
  samples_dist = vegdist(t(unlabeled_proteins_unmelted), method = "bray")
  
  # Caluclate PCoA on distance matrix
  pcoa = pcoa(samples_dist)
  axis_1_variance = round(pcoa$values$Relative_eig[1]*100,digits = 2)
  axis_2_variance = round(pcoa$values$Relative_eig[2]*100, digits = 2)
  pcoa_vectors = data.frame(pcoa$vectors[,1:2])
  vectors_to_plot = pcoa_vectors
  colnames(vectors_to_plot) = c("Axis.1","Axis.2")
  vectors_to_plot$sample = rownames(vectors_to_plot)
  vectors_to_plot$treatment = substr(vectors_to_plot$sample,1,nchar(as.vector(vectors_to_plot$sample))-1)
  
  fig1c_plot = ggplot(data=vectors_to_plot, aes(x=Axis.1, y = Axis.2,fill=treatment, color=treatment)) +
    geom_hline(yintercept = 0, lty="dotted")+
    geom_vline(xintercept = 0, lty="dotted")+
    geom_polygon(size=0.5,alpha=0.2) +
    geom_point(shape=21, size= 2.5) +
    scale_color_manual(name="",guide = guide_legend(reverse = TRUE),values = c("CEx"="gray80","CLy"="gray60","DEx"="gray40","DLy"="gray20")) +
    scale_fill_manual(name="",guide = guide_legend(reverse = TRUE),values = c("CEx"="gray80","CLy"="gray60","DEx"="gray40","DLy"="gray20")) +
    xlab(paste("PCoA 1 (",axis_1_variance,"%)",sep="")) + 
    ylab(paste("PCoA 2 (",axis_2_variance,"%)",sep="")) + 
    xlim(c(-0.3,0.4)) +
    ylim(c(-0.2,0.14)) +
    theme(panel.border = element_rect(colour = "black", fill=NA),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.title = element_text(size=12,color="black"),
          axis.text.x=element_text(size=12,color="black"), 
          axis.text.y=element_text(size=12,color="black"))
  
  ggsave("fig1c.svg", fig1c_plot, device = "svg",scale = 1, width = 3.5, height = 2.5, units = "in", dpi = 600)
  
}

#####################################
############# FIGURE 2a ############# 
#####################################
if(plot_fig2a == "T"){
  
  # Read in data
  all_data = read.csv("all_protein_data.csv", header=T, sep=",", stringsAsFactors = F, comment.char = "")
  all_data$class = ifelse(grepl("Bacteroidetes", all_data$phylum), "Bacteroidetes",as.character(all_data$class))
  
  # Aggregate rows
  sample_by_enrichment_df = aggregate(. ~ class + order + sample, data=all_data[,c(1,7,8,12:ncol(all_data))], FUN=sum)
  sample_by_enrichment_df = sample_by_enrichment_df[!grepl("^t",sample_by_enrichment_df$sample),]
  
  # Aggregate columns
  sample_by_enrichment_df$sum = rowSums(sample_by_enrichment_df[,c(6:ncol(sample_by_enrichment_df))])
  labeled_proteins = sample_by_enrichment_df[,c(1:3,ncol(sample_by_enrichment_df))]
  labeled_proteins$taxon = paste(labeled_proteins$class,labeled_proteins$order,sep="__")
  labeled_proteins$sample = substr(labeled_proteins$sample,1,nchar(as.vector(labeled_proteins$sample))-1)
  labeled_proteins = aggregate(sum ~ taxon + sample, labeled_proteins, FUN=mean)
  
  # Unmelt data to create relative abundance matrix
  labeled_proteins_unmelted = dcast(data = labeled_proteins,formula = sample~taxon,fun.aggregate = sum,value.var = "sum")
  rownames(labeled_proteins_unmelted) = labeled_proteins_unmelted$sample
  labeled_proteins_unmelted$sample = NULL
  labeled_proteins_unmelted = as.data.frame(t(labeled_proteins_unmelted))
  cols = as.character(colnames(labeled_proteins_unmelted))
  rows = as.character(rownames(labeled_proteins_unmelted))
  labeled_proteins_unmelted = labeled_proteins_unmelted %>% mutate_all(as.character)
  labeled_proteins_unmelted = sapply(labeled_proteins_unmelted, as.numeric)
  labeled_proteins_unmelted = sweep(labeled_proteins_unmelted, 2, colSums(labeled_proteins_unmelted), '/')
  rownames(labeled_proteins_unmelted) = rows
  colnames(labeled_proteins_unmelted) = cols
  labeled_proteins_unmelted = as.data.frame(labeled_proteins_unmelted)
  labeled_proteins_unmelted$taxon = rownames(labeled_proteins_unmelted)
  labeled_proteins_unmelted = cSplit(labeled_proteins_unmelted, "taxon", "__")
                      
  keep_taxa = c("Rhodobacterales","Pelagibacterales","Rhodospirillales","SAR116_cluster","Rhizobiales","Unassigned_Alphaproteobacteria", ## Alphas
           "Flavobacteriales","Sphingobacteriales","Bacteroidales","Cytophagales", "Unassigned_Bacteroidetes",  ## Bacts
           "SAR86_cluster", "Alteromonadales","Cellvibrionales","Oceanospirillales", "Unassigned_Gammaproteobacteria", ## Gammas
           "Unassigned_Euryarchaeota",  ## Archaea
           "Actinobacteria",  ## Actinos
           "Burkholderiales","Methylophilales",  ## Betas
           "Unclassified",   ## Unclassified
           "Other")
  
  # Add "Other" category
  labeled_proteins_unmelted_subset = data.frame(labeled_proteins_unmelted[labeled_proteins_unmelted$taxon_2 %in% keep_taxa,])
  labeled_proteins_unmelted_subset = transform(labeled_proteins_unmelted_subset, taxon_1 = as.character(taxon_1),taxon_2 = as.character(taxon_2),CEx = as.numeric(CEx),CLy = as.numeric(CLy),DEx = as.numeric(DEx),DLy = as.numeric(DLy))
  other = data.frame(rbind(1-colSums(labeled_proteins_unmelted_subset[,1:4])),"Other","Other")
  colnames(other)[5:6] = c("taxon_1","taxon_2")
  labeled_proteins_unmelted_subset = rbind(labeled_proteins_unmelted_subset,other)
  labeled_proteins_unmelted_subset$taxon_1 = ifelse(!grepl("Alphaproteobacteria|Gammaproteobacteria|Bacteroidetes", labeled_proteins_unmelted_subset$taxon_1),"Other",labeled_proteins_unmelted_subset$taxon_1)
  labeled_proteins_unmelted_subset$sum = rowSums(labeled_proteins_unmelted_subset[,1:4])
  labeled_proteins_unmelted_subset = labeled_proteins_unmelted_subset[order(-labeled_proteins_unmelted_subset$sum),]
  labeled_proteins_unmelted_subset$sum = NULL

  # Create barplot aes
  barplot_data = split(labeled_proteins_unmelted_subset,f = labeled_proteins_unmelted_subset$taxon_1)
  barplot_data$Alphaproteobacteria$colors = brewer.pal(length(as.character(unique(barplot_data$Alphaproteobacteria$taxon_2))), "Blues")
  barplot_data$Bacteroidetes$colors = brewer.pal(length(as.character(unique(barplot_data$Bacteroidetes$taxon_2))), "Greys")
  barplot_data$Gammaproteobacteria$colors = brewer.pal(length(as.character(unique(barplot_data$Gammaproteobacteria$taxon_2))), "Oranges")
  barplot_data$Other$colors = brewer.pal(length(as.character(unique(barplot_data$Other$taxon_2))), "Purples")

  # Melt data into useable form
  barplot_data$Alphaproteobacteria$taxon_2 = factor(barplot_data$Alphaproteobacteria$taxon_2, levels = barplot_data$Alphaproteobacteria$taxon_2)
  barplot_data$Alphaproteobacteria = melt(barplot_data$Alphaproteobacteria)
  barplot_data$Alphaproteobacteria$variable = factor(barplot_data$Alphaproteobacteria$variable, levels = c("DLy","DEx","CLy","CEx"))
  barplot_data$Bacteroidetes$taxon_2 = factor(barplot_data$Bacteroidetes$taxon_2, levels = barplot_data$Bacteroidetes$taxon_2)
  barplot_data$Bacteroidetes = melt(barplot_data$Bacteroidetes)
  barplot_data$Bacteroidetes$variable = factor(barplot_data$Bacteroidetes$variable, levels = c("DLy","DEx","CLy","CEx"))
  barplot_data$Gammaproteobacteria$taxon_2 = factor(barplot_data$Gammaproteobacteria$taxon_2, levels = barplot_data$Gammaproteobacteria$taxon_2)
  barplot_data$Gammaproteobacteria = melt(barplot_data$Gammaproteobacteria)
  barplot_data$Gammaproteobacteria$variable = factor(barplot_data$Gammaproteobacteria$variable, levels = c("DLy","DEx","CLy","CEx"))
  barplot_data$Other$taxon_2 = factor(barplot_data$Other$taxon_2, levels = barplot_data$Other$taxon_2)
  barplot_data$Other = melt(barplot_data$Other)
  barplot_data$Other$variable = factor(barplot_data$Other$variable, levels = c("DLy","DEx","CLy","CEx"))
  
  barplot1 = ggplot(data=barplot_data$Alphaproteobacteria, aes(x=variable, y=value, group=rev(taxon_2), fill=taxon_2)) +
    geom_bar(show.legend = T,stat="identity", position = position_stack(reverse=TRUE),colour = "gray40", width=0.9) +
    scale_y_continuous(expand = c(0.05,0), limits=c(0,0.8)) +
    scale_fill_manual(guide = "legend", values = barplot_data$Alphaproteobacteria$colors) +
    theme_minimal() +
    xlab("") +
    ylab("") +
    guides(fill=guide_legend(title="New Legend Title")) +
    theme(axis.text.y = element_text(size=10, colour="black"),
          axis.text.x = element_text(size = 10, colour = "black"), 
          axis.title.x=element_blank(),
          legend.text = element_text(size=10),
          legend.key.size = unit(0.6,"cm"),
          legend.title = element_blank())
  
  # Define other barplots based on the original
  barplot2 = barplot1 %+% barplot_data$Bacteroidetes + xlab("Relative abundance in metaproteome") + scale_fill_manual(guide = "legend", values = barplot_data$Bacteroidetes$colors)
  barplot3 = barplot1 %+% barplot_data$Gammaproteobacteria + scale_fill_manual(guide = "legend", values = barplot_data$Gammaproteobacteria$colors)
  barplot4 = barplot1 %+% barplot_data$Other + scale_fill_manual(guide = "legend", values = barplot_data$Other$colors)
  
  # Extract legends
  barplot1_legend = get_legend(barplot1)
  barplot2_legend = get_legend(barplot2)
  barplot3_legend = get_legend(barplot3)
  barplot4_legend = get_legend(barplot4)
  
  # Remove legends for first facet
  barplot1 = barplot1 + theme(legend.position = "none") 
  barplot2 = barplot2 + theme(legend.position = "none") 
  barplot3 = barplot3 + theme(legend.position = "none") 
  barplot4 = barplot4 + theme(legend.position = "none") 
  
  # Create barplot facet and legend facet
  all_barplots = grid.arrange(barplot1, barplot2, barplot3, ncol = 1)
  all_barplot_legends = grid.arrange(barplot1_legend, barplot2_legend, barplot3_legend, ncol = 1)

  # Stick together bars and legends
  fig2a_plot = grid.arrange(all_barplots,all_barplot_legends,ncol=2)
  
  # Print
  ggsave("fig2a.svg", fig2a_plot, device = "svg",scale = 1, width = 5, height = 5, units = "in", dpi = 600)
  
}

#####################################
############# FIGURE 2b ############# 
#####################################
if(plot_fig2b == "T"){
  
  # Read in data
  all_data = read.csv("all_protein_data.csv", header=T, sep=",", stringsAsFactors = F, comment.char = "")
  all_data$class = ifelse(grepl("Bacteroidetes", all_data$phylum), "Bacteroidetes",as.character(all_data$class))
  keep_taxa = c("Alphaproteobacteria","Bacteroidetes","Gammaproteobacteria")
  all_data = all_data[all_data$class %in% keep_taxa,]
  
  # Remove extraneous information
  sample_by_class_by_enrichment_df = aggregate(. ~ sample + class, data=all_data[,c(1,7,12:ncol(all_data))], FUN=sum)
  sample_by_class_by_enrichment_df = sample_by_class_by_enrichment_df[!grepl("^t",sample_by_class_by_enrichment_df$sample),]
  
  # Define labeled (>= 2%) and unlabeled (< 2%) fractions 
  samples = sample_by_class_by_enrichment_df$sample
  classes = sample_by_class_by_enrichment_df$class
  all_proteins = sample_by_class_by_enrichment_df[,3:ncol(sample_by_class_by_enrichment_df)]
  unlabeled_proteins = sample_by_class_by_enrichment_df[,3:4]
  labeled_proteins = sample_by_class_by_enrichment_df[,5:ncol(sample_by_class_by_enrichment_df)]
  
  # Calculate ENRV
  enrv_multiplier = as.numeric(2:100)
  enrv = data.frame(mapply(`*`,labeled_proteins,enrv_multiplier))
  enrv = rowSums(enrv)/rowSums(labeled_proteins)
  
  # Calculate LF
  lf = (rowSums(labeled_proteins)/rowSums(all_proteins))*100
  
  # Prepare data for plotting
  enrv_lf_summary = data.frame(samples,classes,enrv,lf)
  colnames(enrv_lf_summary) = c("sample","class","enrv","lf")
  
  # Calculate histogram for plotting
  master_hist = as.data.frame(t(apply(sample_by_class_by_enrichment_df[,5:ncol(sample_by_class_by_enrichment_df)], 1, tapply, gl(33, 3, ncol(sample_by_class_by_enrichment_df[,5:ncol(sample_by_class_by_enrichment_df)])), mean)))
  master_hist$sample = samples
  master_hist$class  = classes
  master_hist$sample = substr(master_hist$sample,1,nchar(as.vector(master_hist$sample))-1)
  master_hist_mean = aggregate(. ~ sample + class, data=master_hist, FUN=mean)
  master_hist_mean_melted = melt(master_hist_mean)
  master_hist_mean_melted$value = log10(master_hist_mean_melted$value+1)
  
  # Prepare for plotting
  master_hist_mean_melted_split = split(master_hist_mean_melted,f = master_hist_mean_melted$class)
  
  hist1 = ggplot(data=master_hist_mean_melted_split$Alphaproteobacteria, aes(x=variable, y=value, linetype=sample,color=sample, group = sample)) +
    geom_line(size=0.6, show.legend=T) +
    scale_color_manual(name="",guide= "legend",values = c("CEx"="gray80","CLy"="gray60","DEx"="gray40","DLy"="gray20")) +
    scale_linetype_manual(name="",guide= "legend",values = c("CEx"="twodash","CLy"="solid","DEx"="twodash","DLy"="solid")) +
    theme_minimal()+
    scale_x_discrete(labels = c("2",rep("",15),"50",rep("",15),"100")) +
    scale_y_continuous(expand=c(0,0),limits= c(0,2.7), labels=c(10^0,10^1,10^2,"500"), breaks=c(0,1,2,2.7)) +    
    ggtitle("Alphaproteobacteria") +
    xlab("") +
    ylab("") +
    theme(axis.text.x = element_text(size=10),
          axis.text.y = element_text(size=10),
          legend.title = element_blank(),
          strip.background = element_blank(),
          strip.text.x = element_blank(),
          plot.title = element_text(hjust = 0.5, size=13,face="bold"),
          axis.line.y = element_line(color = "black"),
          axis.line.x = element_line(color = "black"),
          axis.ticks = element_line(color = "black"),
          axis.text = element_text(color = "black"),
          axis.title = element_text(color = "black"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
  
  hist2 = hist1 %+% master_hist_mean_melted_split$Bacteroidetes
  hist2 = hist2 + ggtitle("Bacteroidetes") + ylab("NBSC")
  
  hist3 = hist1 %+% master_hist_mean_melted_split$Gammaproteobacteria
  hist3 = hist3 + ggtitle("Gammaproteobacteria") + xlab("ENRV (3% bins)")
  
  # Extract legend
  hist_legend = get_legend(hist1)
  hist1 = hist1 + theme(legend.position = "none")
  hist2 = hist2 + theme(legend.position = "none")
  hist3 = hist3 + theme(legend.position = "none")
  
  # Combine hists and legend
  fig2_plot = grid.arrange(hist1, hist2, hist3, ncol=1)
  fig2_plot = grid.arrange(fig2_plot, hist_legend, ncol=2)
  
  ggsave("fig2b.svg", fig2_plot, device = "svg",scale = 1, width = 5, height = 5, units = "in", dpi = 600)
  
}

#####################################
############# FIGURE 3a ############# 
#####################################
if(plot_fig3a == "T"){
  
  # Read in data
  all_data = read.csv("all_protein_data.csv", header=T, sep=",", stringsAsFactors = F, comment.char = "")
  all_data$class = ifelse(grepl("Bacteroidetes", all_data$phylum), "Bacteroidetes",as.character(all_data$class))
  
  # Define plotting variables based on abundant enzymes/pathways/functions
  all_data$cog_desc = ifelse(grepl("^Ribosomal_protein",all_data$cog_desc),"Ribosomal_proteins",as.character(all_data$cog_desc))
  all_data$cog_desc = ifelse(grepl("polymerase",all_data$cog_desc),"DNA/RNA_polymerases",as.character(all_data$cog_desc))
  all_data$cog_desc = ifelse(grepl("^F0F1",all_data$cog_desc),"ATPases",as.character(all_data$cog_desc))
  all_data$cog_desc = ifelse(grepl("vacuolar-type",all_data$cog_desc),"ATPases",as.character(all_data$cog_desc))
  all_data$cog_desc = ifelse(grepl("^MoxR",all_data$cog_desc),"ATPases",as.character(all_data$cog_desc))
  all_data$cog_desc = ifelse(grepl("ranscription",all_data$cog_desc),"Transcription/Translation_factors",as.character(all_data$cog_desc))
  all_data$cog_desc = ifelse(grepl("ranslation",all_data$cog_desc),"Transcription/Translation_factors",as.character(all_data$cog_desc))
  all_data$cog_desc = ifelse(grepl("Polyribonucleotide",all_data$cog_desc),"Transcription/Translation_factors",as.character(all_data$cog_desc))
  all_data$cog_desc = ifelse(grepl("Biopolymer_transport_protein",all_data$cog_desc),"Biopolymer_transport",as.character(all_data$cog_desc))
  all_data$cog_desc = ifelse(grepl("Flagellar_motor_protein",all_data$cog_desc),"Biopolymer_transport",as.character(all_data$cog_desc))
  all_data$cog_desc = ifelse(grepl("Periplasmic_protein_involved_in_polysaccharide_export",all_data$cog_desc),"Biopolymer_transport",as.character(all_data$cog_desc))
  all_data$cog_desc = ifelse(grepl("Poly(3-hydroxybutyrate)_depolymerase",all_data$cog_desc),"Biopolymer_transport",as.character(all_data$cog_desc))
  all_data$cog_desc = ifelse(grepl("permidine",all_data$cog_desc),"Polyamine_transport",as.character(all_data$cog_desc))
  all_data$cog_desc = ifelse(grepl("ABC-type",all_data$cog_desc),"ABC_transport",as.character(all_data$cog_desc))
  all_data$cog_desc = ifelse(grepl("ABC_trans",all_data$cog_desc),"ABC_transport",as.character(all_data$cog_desc))
  all_data$cog_desc = ifelse(grepl("TRAP-type",all_data$cog_desc),"TRAP_transport",as.character(all_data$cog_desc))
  all_data$cog_desc = ifelse(grepl("Outer_membrane_receptor_proteins_",all_data$cog_desc),"TonB_transport",as.character(all_data$cog_desc))
  all_data$cog_desc = ifelse(grepl("TonB",all_data$cog_desc),"TonB_transport",as.character(all_data$cog_desc))
  all_data$cog_desc = ifelse(grepl("Outer_membrane_cobalamin_receptor_protein",all_data$cog_desc),"TonB_transport",as.character(all_data$cog_desc))
  all_data$cog_desc = ifelse(grepl("Outer_membrane_receptor_for_ferrienterochelin_and_colicins",all_data$cog_desc),"TonB_transport",as.character(all_data$cog_desc))
  all_data$cog_desc = ifelse(grepl("Riboso",all_data$cog_desc),"Ribosomal_proteins",as.character(all_data$cog_desc))
  all_data$cog_desc = gsub('_system_periplasmic_component', '', as.character(all_data$cog_desc))
  all_data$cog_desc = gsub('_systems_periplasmic_component', '', as.character(all_data$cog_desc))
  all_data$cog_desc = gsub('/domain', '', as.character(all_data$cog_desc))
  
  # Aggregate rows
  sample_by_enrichment_df = aggregate(. ~ cog_desc + class + sample, data=all_data[,c(1,4,7,12:ncol(all_data))], FUN=sum)
  sample_by_enrichment_df = sample_by_enrichment_df[grepl("^t",sample_by_enrichment_df$sample),]
  keep_classes = c("Alphaproteobacteria","Bacteroidetes","Gammaproteobacteria")
  keep_functions  = c("Nucleoside_diphosphate_kinase","TonB_transport","ABC_transport","TRAP_transport","Ribosomal_proteins", "Transcription/Translation_factors", "ATPases", "DNA/RNA_polymerases", "Polyamine_transport","Biopolymer_transport")
  sample_by_enrichment_df = sample_by_enrichment_df[sample_by_enrichment_df$class %in% keep_classes,]
  sample_by_enrichment_df = sample_by_enrichment_df[sample_by_enrichment_df$cog_desc %in% keep_functions,]
  
  # Aggregate columns
  sample_by_enrichment_df$sum = rowSums(sample_by_enrichment_df[,4:5])
  unlabeled_proteins = sample_by_enrichment_df[,c(1:3,ncol(sample_by_enrichment_df))]
  unlabeled_proteins$taxon__function = paste(unlabeled_proteins$class,unlabeled_proteins$cog_desc,sep="__")
  unlabeled_proteins$sample = substr(unlabeled_proteins$sample,1,nchar(as.vector(unlabeled_proteins$sample))-1)
  unlabeled_proteins = aggregate(sum ~ taxon__function + sample, unlabeled_proteins, FUN=mean)
  
  # Unmelt data to create relative abundance matrix
  unlabeled_proteins_unmelted = dcast(data = unlabeled_proteins,formula = sample~taxon__function,fun.aggregate = sum,value.var = "sum")
  rownames(unlabeled_proteins_unmelted) = unlabeled_proteins_unmelted$sample
  unlabeled_proteins_unmelted$sample = NULL
  unlabeled_proteins_unmelted = as.data.frame(t(unlabeled_proteins_unmelted))
  cols = as.character(colnames(unlabeled_proteins_unmelted))
  rows = as.character(rownames(unlabeled_proteins_unmelted))
  unlabeled_proteins_unmelted = unlabeled_proteins_unmelted %>% mutate_all(as.character)
  unlabeled_proteins_unmelted = sapply(unlabeled_proteins_unmelted, as.numeric)
  rownames(unlabeled_proteins_unmelted) = rows
  colnames(unlabeled_proteins_unmelted) = cols
  unlabeled_proteins_unmelted = as.data.frame(unlabeled_proteins_unmelted)
  unlabeled_proteins_unmelted$taxon__function = rownames(unlabeled_proteins_unmelted)
  unlabeled_proteins_unmelted = cSplit(unlabeled_proteins_unmelted, "taxon__function", "__")
  unlabeled_proteins_remelted = melt(unlabeled_proteins_unmelted)
  
  # Split into samples for plotting 
  ternary_input_data = split(unlabeled_proteins_remelted,f = unlabeled_proteins_remelted$variable)
  
  # Define T0 values and legend
  ternary_input_data$t0 = dcast(data = ternary_input_data$t0, formula = taxon__function_2 ~ taxon__function_1, fun.aggregate = sum,value.var = "value")
  t0_function_sums = rowSums(ternary_input_data$t0[,2:ncol(ternary_input_data$t0)])
  number_to_function_lookup = data.frame(rownames(ternary_input_data$t0),ternary_input_data$t0$taxon__function_2)
  colnames(number_to_function_lookup) = c("Label","Function")
  number_to_function_lookup$Function = gsub("_"," ",number_to_function_lookup$Function)
  rownames(number_to_function_lookup) = number_to_function_lookup$Label
  number_to_function_lookup$Label = NULL
  ternary_input_data$t0$taxon__function_2 = NULL
  ternary_input_data$t0 = t(ternary_input_data$t0)
  ternary_input_data$t0 = as.data.frame(t(sweep(ternary_input_data$t0, 2, colSums(ternary_input_data$t0), '/')))
  
  # Define T1 values
  ternary_input_data$t1 = dcast(data = ternary_input_data$t1, formula = taxon__function_2 ~ taxon__function_1, fun.aggregate = sum,value.var = "value")
  t1_function_sums = rowSums(ternary_input_data$t1[,2:ncol(ternary_input_data$t1)])
  ternary_input_data$t1$taxon__function_2 = NULL
  ternary_input_data$t1 = t(ternary_input_data$t1)
  ternary_input_data$t1 = as.data.frame(t(sweep(ternary_input_data$t1, 2, colSums(ternary_input_data$t1), '/')))
  
  # Set up file for plotting
  svg(filename = "fig3a.svg",width = 2.5, height = 6)
  par(mfrow=c(2, 1), mar=rep(0, 4))
  makeTransparent = function(someColor, alpha=100){ newColor<-col2rgb(someColor); apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],blue=curcoldata[3],alpha=alpha, maxColorValue=255)})}
  
  # T0 plot
  TernaryPlot(main = "T0", alab="Alphaproteobacteria", blab="Bacteroidetes", clab="Gammaproteobacteria",lab.col=c('#2171B5', 'gray50', 'darkorange'),lab.cex=1, grid.minor.lines=0,grid.lty='solid', padding=0.14,grid.col='gray93', axis.col = "gray50",ticks.col="gray50")
  AddToTernary(points, ternary_input_data$t0, pch=21, col="black",bg=makeTransparent('gray50',alpha = 100),cex=log10(t0_function_sums+1))
  AddToTernary(text, ternary_input_data$t0, as.character(rownames(ternary_input_data$t0)), cex=0.6, font=2)
  
  # T1 plot
  TernaryPlot(main = "T1", alab="Alphaproteobacteria", blab="Bacteroidetes", clab="Gammaproteobacteria",lab.col=c('#2171B5', 'gray50', 'darkorange'),lab.cex=1, grid.minor.lines=0,grid.lty='solid', padding=0.14,grid.col='gray93', axis.col = "gray50",ticks.col="gray50")
  AddToTernary(points, ternary_input_data$t1, pch=21, col="black",bg=makeTransparent('gray50',alpha = 100),cex=log10(t1_function_sums+1))
  AddToTernary(text, ternary_input_data$t1, as.character(rownames(ternary_input_data$t1)), cex=0.6, font=2)
  
  dev.off()
  
  # Plot legend
  svg(filename = "fig3c.svg",width = 3, height = 3)
  legend = tableGrob(number_to_function_lookup)
  grid.arrange(legend)
  dev.off()
}

#####################################
############# FIGURE 3b ############# 
#####################################
if(plot_fig3b == "T"){

  # Read in data
  all_data = read.csv("all_protein_data.csv", header=T, sep=",", stringsAsFactors = F, comment.char = "")
  all_data$class = ifelse(grepl("Bacteroidetes", all_data$phylum), "Bacteroidetes",as.character(all_data$class))
  
  # Define plotting variables based on abundant enzymes/pathways/functions
  all_data$cog_desc = ifelse(grepl("^COG3181",all_data$cog),"Tripartite_tricarboxylate_transport",as.character(all_data$cog_desc))
  all_data$cog_desc = ifelse(grepl("^Ribosomal_protein",all_data$cog_desc),"Ribosomal_proteins",as.character(all_data$cog_desc))
  all_data$cog_desc = ifelse(grepl("polymerase",all_data$cog_desc),"DNA/RNA_polymerases",as.character(all_data$cog_desc))
  all_data$cog_desc = ifelse(grepl("^F0F1",all_data$cog_desc),"ATPases",as.character(all_data$cog_desc))
  all_data$cog_desc = ifelse(grepl("vacuolar-type",all_data$cog_desc),"ATPases",as.character(all_data$cog_desc))
  all_data$cog_desc = ifelse(grepl("^MoxR",all_data$cog_desc),"ATPases",as.character(all_data$cog_desc))
  all_data$cog_desc = ifelse(grepl("ranscription",all_data$cog_desc),"Transcription/Translation_factors",as.character(all_data$cog_desc))
  all_data$cog_desc = ifelse(grepl("ranslation",all_data$cog_desc),"Transcription/Translation_factors",as.character(all_data$cog_desc))
  all_data$cog_desc = ifelse(grepl("Polyribonucleotide",all_data$cog_desc),"Transcription/Translation_factors",as.character(all_data$cog_desc))
  all_data$cog_desc = ifelse(grepl("Biopolymer_transport_protein",all_data$cog_desc),"Biopolymer_transport",as.character(all_data$cog_desc))
  all_data$cog_desc = ifelse(grepl("Flagellar_motor_protein",all_data$cog_desc),"Biopolymer_transport",as.character(all_data$cog_desc))
  all_data$cog_desc = ifelse(grepl("Periplasmic_protein_involved_in_polysaccharide_export",all_data$cog_desc),"Biopolymer_transport",as.character(all_data$cog_desc))
  all_data$cog_desc = ifelse(grepl("Poly(3-hydroxybutyrate)_depolymerase",all_data$cog_desc),"Biopolymer_transport",as.character(all_data$cog_desc))
  all_data$cog_desc = ifelse(grepl("permidine",all_data$cog_desc),"Polyamine_transport",as.character(all_data$cog_desc))
  all_data$cog_desc = ifelse(grepl("ABC-type",all_data$cog_desc),"ABC_transport",as.character(all_data$cog_desc))
  all_data$cog_desc = ifelse(grepl("ABC_trans",all_data$cog_desc),"ABC_transport",as.character(all_data$cog_desc))
  all_data$cog_desc = ifelse(grepl("TRAP-type",all_data$cog_desc),"TRAP_transport",as.character(all_data$cog_desc))
  all_data$cog_desc = ifelse(grepl("Outer_membrane_receptor_proteins_",all_data$cog_desc),"TonB_transport",as.character(all_data$cog_desc))
  all_data$cog_desc = ifelse(grepl("TonB",all_data$cog_desc),"TonB_transport",as.character(all_data$cog_desc))
  all_data$cog_desc = ifelse(grepl("Outer_membrane_cobalamin_receptor_protein",all_data$cog_desc),"TonB_transport",as.character(all_data$cog_desc))
  all_data$cog_desc = ifelse(grepl("Outer_membrane_receptor_for_ferrienterochelin_and_colicins",all_data$cog_desc),"TonB_transport",as.character(all_data$cog_desc))
  all_data$cog_desc = ifelse(grepl("Riboso",all_data$cog_desc),"Ribosomal_proteins",as.character(all_data$cog_desc))
  all_data$cog_desc = gsub('_system_periplasmic_component', '', as.character(all_data$cog_desc))
  all_data$cog_desc = gsub('_systems_periplasmic_component', '', as.character(all_data$cog_desc))
  all_data$cog_desc = gsub('/domain', '', as.character(all_data$cog_desc))

  # Aggregate rows
  sample_by_enrichment_df = aggregate(. ~ cog_desc + class + sample, data=all_data[,c(1,4,7,12:ncol(all_data))], FUN=sum)
  sample_by_enrichment_df = sample_by_enrichment_df[!grepl("^t",sample_by_enrichment_df$sample),]
  keep_classes = c("Alphaproteobacteria","Bacteroidetes","Gammaproteobacteria")
  keep_functions  = c("Nucleoside_diphosphate_kinase","TonB_transport","ABC_transport","TRAP_transport","Ribosomal_proteins", "Transcription/Translation_factors", "ATPases", "DNA/RNA_polymerases", "Polyamine_transport","Biopolymer_transport")
  sample_by_enrichment_df = sample_by_enrichment_df[sample_by_enrichment_df$class %in% keep_classes,]
  sample_by_enrichment_df = sample_by_enrichment_df[sample_by_enrichment_df$cog_desc %in% keep_functions,]
  
  # Aggregate columns
  sample_by_enrichment_df$sum = rowSums(sample_by_enrichment_df[,c(6:ncol(sample_by_enrichment_df))])
  labeled_proteins = sample_by_enrichment_df[,c(1:3,ncol(sample_by_enrichment_df))]
  labeled_proteins$taxon__function = paste(labeled_proteins$class,labeled_proteins$cog_desc,sep="__")
  labeled_proteins$sample = substr(labeled_proteins$sample,1,nchar(as.vector(labeled_proteins$sample))-1)
  labeled_proteins = aggregate(sum ~ taxon__function + sample, labeled_proteins, FUN=mean)
  
  # Unmelt data to create relative abundance matrix
  labeled_proteins_unmelted = dcast(data = labeled_proteins,formula = sample~taxon__function,fun.aggregate = sum,value.var = "sum")
  rownames(labeled_proteins_unmelted) = labeled_proteins_unmelted$sample
  labeled_proteins_unmelted$sample = NULL
  labeled_proteins_unmelted = as.data.frame(t(labeled_proteins_unmelted))
  cols = as.character(colnames(labeled_proteins_unmelted))
  rows = as.character(rownames(labeled_proteins_unmelted))
  labeled_proteins_unmelted = labeled_proteins_unmelted %>% mutate_all(as.character)
  labeled_proteins_unmelted = sapply(labeled_proteins_unmelted, as.numeric)
  rownames(labeled_proteins_unmelted) = rows
  colnames(labeled_proteins_unmelted) = cols
  labeled_proteins_unmelted = as.data.frame(labeled_proteins_unmelted)
  labeled_proteins_unmelted$taxon__function = rownames(labeled_proteins_unmelted)
  labeled_proteins_unmelted = cSplit(labeled_proteins_unmelted, "taxon__function", "__")
  labeled_proteins_remelted = melt(labeled_proteins_unmelted)
  
  # Split into samples for plotting
  ternary_input_data = split(labeled_proteins_remelted,f = labeled_proteins_remelted$variable)
  
  # Define CEx values and legend
  ternary_input_data$CEx = dcast(data = ternary_input_data$CEx, formula = taxon__function_2 ~ taxon__function_1, fun.aggregate = sum,value.var = "value")
  CEx_function_sums = rowSums(ternary_input_data$CEx[,2:ncol(ternary_input_data$CEx)])
  number_to_function_lookup = data.frame(rownames(ternary_input_data$CEx),ternary_input_data$CEx$taxon__function_2)
  colnames(number_to_function_lookup) = c("Label","Function")
  number_to_function_lookup$Function = gsub("_"," ",number_to_function_lookup$Function)
  rownames(number_to_function_lookup) = number_to_function_lookup$Label
  number_to_function_lookup$Label = NULL
  ternary_input_data$CEx$taxon__function_2 = NULL
  ternary_input_data$CEx = t(ternary_input_data$CEx)
  ternary_input_data$CEx = as.data.frame(t(sweep(ternary_input_data$CEx, 2, colSums(ternary_input_data$CEx), '/')))

  # Define CLy values
  ternary_input_data$CLy = dcast(data = ternary_input_data$CLy, formula = taxon__function_2 ~ taxon__function_1, fun.aggregate = sum,value.var = "value")
  CLy_function_sums = rowSums(ternary_input_data$CLy[,2:ncol(ternary_input_data$CLy)])
  ternary_input_data$CLy$taxon__function_2 = NULL
  ternary_input_data$CLy = t(ternary_input_data$CLy)
  ternary_input_data$CLy = as.data.frame(t(sweep(ternary_input_data$CLy, 2, colSums(ternary_input_data$CLy), '/')))
  
  # Define DEx values
  ternary_input_data$DEx = dcast(data = ternary_input_data$DEx, formula = taxon__function_2 ~ taxon__function_1, fun.aggregate = sum,value.var = "value")
  DEx_function_sums = rowSums(ternary_input_data$DEx[,2:ncol(ternary_input_data$DEx)])
  ternary_input_data$DEx$taxon__function_2 = NULL
  ternary_input_data$DEx = t(ternary_input_data$DEx)
  ternary_input_data$DEx = as.data.frame(t(sweep(ternary_input_data$DEx, 2, colSums(ternary_input_data$DEx), '/')))

  # Define DLy values
  ternary_input_data$DLy = dcast(data = ternary_input_data$DLy, formula = taxon__function_2 ~ taxon__function_1, fun.aggregate = sum,value.var = "value")
  DLy_function_sums = rowSums(ternary_input_data$DLy[,2:ncol(ternary_input_data$DLy)])
  ternary_input_data$DLy$taxon__function_2 = NULL
  ternary_input_data$DLy = t(ternary_input_data$DLy)
  ternary_input_data$DLy = as.data.frame(t(sweep(ternary_input_data$DLy, 2, colSums(ternary_input_data$DLy), '/')))
  
  # Set up plot file
  svg(filename = "fig3b.svg",width = 5, height = 6)
  par(mfrow=c(2, 2), mar=rep(0, 4))
  
  # Plot CEx
  TernaryPlot(main = "CEx", alab="Alphaproteobacteria", blab="Bacteroidetes", clab="Gammaproteobacteria",lab.col=c('#2171B5', 'gray50', 'darkorange'),lab.cex=1.2, grid.minor.lines=0,grid.lty='solid', padding=0.14,grid.col='gray93', axis.col = "gray50",ticks.col="gray50")
  AddToTernary(points, ternary_input_data$CEx, pch=21, col="black",bg=makeTransparent('gray50',alpha = 100),cex=log10(CEx_function_sums+1))
  AddToTernary(text, ternary_input_data$CEx, as.character(rownames(ternary_input_data$CEx)), cex=0.6, font=2)
  
  # Plot CLy
  TernaryPlot(main = "CLy", alab="Alphaproteobacteria", blab="Bacteroidetes", clab="Gammaproteobacteria",lab.col=c('#2171B5', 'gray50', 'darkorange'),lab.cex=1.2, grid.minor.lines=0,grid.lty='solid', padding=0.14,grid.col='gray93', axis.col = "gray50",ticks.col="gray50")
  AddToTernary(points, ternary_input_data$CLy, pch=21, col="black",bg=makeTransparent('gray50',alpha = 100),cex=log10(CLy_function_sums+1))
  AddToTernary(text, ternary_input_data$CLy, as.character(rownames(ternary_input_data$CLy)), cex=0.6, font=2)
  
  # Plot DEx
  TernaryPlot(main = "DEx", alab="Alphaproteobacteria", blab="Bacteroidetes", clab="Gammaproteobacteria",lab.col=c('#2171B5', 'gray50', 'darkorange'),lab.cex=1.2, grid.minor.lines=0,grid.lty='solid', padding=0.14,grid.col='gray93', axis.col = "gray50",ticks.col="gray50")
  AddToTernary(points, ternary_input_data$DEx, pch=21, col="black",bg=makeTransparent('gray50',alpha = 100),cex=log10(DEx_function_sums+1))
  AddToTernary(text, ternary_input_data$DEx, as.character(rownames(ternary_input_data$DEx)), cex=0.6, font=2)
  
  # Plot DLy
  TernaryPlot(main = "DLy", alab="Alphaproteobacteria", blab="Bacteroidetes", clab="Gammaproteobacteria",lab.col=c('#2171B5', 'gray50', 'darkorange'),lab.cex=1.2, grid.minor.lines=0,grid.lty='solid', padding=0.14,grid.col='gray93', axis.col = "gray50",ticks.col="gray50")
  AddToTernary(points, ternary_input_data$DLy, pch=21, col="black",bg=makeTransparent('gray50',alpha = 100),cex=log10(DLy_function_sums+1))
  AddToTernary(text, ternary_input_data$DLy, as.character(rownames(ternary_input_data$DLy)), cex=0.6, font=2)
  
  dev.off()
  
  # Plot legend
  svg(filename = "fig3c.svg",width = 3, height = 3)
  legend = tableGrob(number_to_function_lookup)
  grid.arrange(legend)
  dev.off()
  
}

#####################################
############# FIGURE 4a ############# 
#####################################
if(plot_fig4a == "T"){
  
  # Read in data
  all_data = read.csv("all_protein_data.csv", header=T, sep=",", stringsAsFactors = F, comment.char = "")
  all_data$class = ifelse(grepl("Bacteroidetes", all_data$phylum), "Bacteroidetes",as.character(all_data$class))
  keep_taxa = c("Alphaproteobacteria","Bacteroidetes","Gammaproteobacteria")
  all_data = all_data[all_data$class %in% keep_taxa,]
  
  # Remove extraneous information
  sample_by_class_by_enrichment_df = aggregate(. ~ sample + class, data=all_data[,c(1,7,12:ncol(all_data))], FUN=sum)
  sample_by_class_by_enrichment_df = sample_by_class_by_enrichment_df[!grepl("^t",sample_by_class_by_enrichment_df$sample),]
  
  # Define labeled (>= 2%) and unlabeled (< 2%) fractions 
  samples = sample_by_class_by_enrichment_df$sample
  classes = sample_by_class_by_enrichment_df$class
  all_proteins = sample_by_class_by_enrichment_df[,3:ncol(sample_by_class_by_enrichment_df)]
  unlabeled_proteins = sample_by_class_by_enrichment_df[,3:4]
  labeled_proteins = sample_by_class_by_enrichment_df[,5:ncol(sample_by_class_by_enrichment_df)]
  
  # Calculate ENRV
  enrv_multiplier = as.numeric(2:100)
  enrv = data.frame(mapply(`*`,labeled_proteins,enrv_multiplier))
  enrv = rowSums(enrv)/rowSums(labeled_proteins)
  
  # Calculate LF
  lf = (rowSums(labeled_proteins)/rowSums(all_proteins))*100
  
  # Prepare data for plotting
  enrv_lf_summary = data.frame(samples,classes,enrv*lf)
  colnames(enrv_lf_summary) = c("sample","class","enrv_times_lf")
  enrv_lf_summary$sample = substr(enrv_lf_summary$sample,1,nchar(as.vector(enrv_lf_summary$sample))-1)
  
  # Calculate means and standard errors
  means = aggregate(enrv_times_lf ~ sample + class, enrv_lf_summary, mean)
  ses = aggregate(enrv_times_lf ~ sample + class, enrv_lf_summary, sd)
  plotting_data = cbind(means,ses$enrv_times_lf/sqrt(3))
  colnames(plotting_data) = c("sample","class","means","ses")
  plotting_data$means = sqrt(plotting_data$means)
  plotting_data$ses = sqrt(plotting_data$ses)
  plotting_data$sample = factor(plotting_data$sample, levels = c("DLy","DEx","CLy","CEx"))
  
  # Plot
  fig4a_plot = ggplot(plotting_data, aes(x=sample,y=means,group=class,fill=class)) +
    geom_errorbar(position=position_dodge(0.4),aes(ymin=means-ses,ymax=means+ses),color="black",width = 0.2) +
    geom_point(position=position_dodge(0.4),shape=21,size=2.5,color="black",stroke=0.5) +
    scale_fill_manual(values = c("Alphaproteobacteria"='#2171B5', "Bacteroidetes" ='gray50', "Gammaproteobacteria" ='darkorange')) +
    theme_minimal() +
    xlab("") +
    ylab("Assimilation sqrt(ENRV*LF)") +
    theme(panel.border = element_rect(colour = "black", fill=NA),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.title = element_text(size=12,color="black"),
          axis.text.x=element_text(size=12,color="black"), 
          axis.text.y=element_text(size=12,color="black"))
  
  ggsave("fig4a.svg", fig4a_plot, device = "svg",scale = 1, width = 4.5, height = 2.5, units = "in", dpi = 600)
  
  
}
  
#####################################
############# FIGURE 4b ############# 
#####################################
if(plot_fig4b == "T"){
  
  # Read in data
  bin_growth_rates = read.table("growth_rate_summary.tsv", header=T, sep="\t", stringsAsFactors = F)
  bin_growth_rates$class = ifelse(grepl("Flavo", bin_growth_rates$class), "Bacteroidetes",as.character(bin_growth_rates$class))
  
  # Calculate means at the class level
  bin_growth_rates_agg = aggregate(. ~ class + sample, bin_growth_rates[,c(1,4:ncol(bin_growth_rates))],mean)

  # Filter just top taxa
  keep_classes = c("Alphaproteobacteria","Bacteroidetes","Gammaproteobacteria")
  bin_growth_rates_agg = bin_growth_rates_agg[bin_growth_rates_agg$class %in% keep_classes,]
  
  # Calculate means and standard errors
  means = rowMeans(bin_growth_rates_agg[,3:5])
  ses = rowSds(bin_growth_rates_agg[,3:5])/sqrt(3)
  plotting_data = cbind(bin_growth_rates_agg[,1:2],means,ses)
  plotting_data = plotting_data[!grepl("^T",plotting_data$sample),]
  plotting_data$sample = factor(plotting_data$sample, levels = c("DLy","DEx","CLy","CEx"))
  
  # Plot
  fig4b_plot = ggplot(plotting_data, aes(x=sample,y=means,group=class,fill=class)) +
    geom_errorbar(position=position_dodge(0.4),aes(ymin=means-ses,ymax=means+ses),color="black",width = 0.2) +
    geom_point(position=position_dodge(0.4),shape=21,size=2.5,color="black",stroke=0.5) +
    scale_fill_manual(values = c("Alphaproteobacteria"='#2171B5', "Bacteroidetes" ='gray50', "Gammaproteobacteria" ='darkorange')) +
    theme_minimal() +
    xlab("") +
    ylab("Estimated Growth Rate") +
    theme(panel.border = element_rect(colour = "black", fill=NA),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.title = element_text(size=12,color="black"),
          axis.text.x=element_text(size=12,color="black"), 
          axis.text.y=element_text(size=12,color="black"))
    
  ggsave("fig4b.svg", fig4b_plot, device = "svg",scale = 1, width = 4.5, height = 2.5, units = "in", dpi = 600)
  
  
}

#####################################
############# FIGURE 4c ############# 
#####################################
if(plot_fig4c == "T"){
  
  # Read in data
  c_to_n_ratio = read.table("c_to_n_ratio_summary.tsv", header=T, sep="\t", stringsAsFactors = F)
  c_to_n_ratio$class = ifelse(grepl("Bacteroidetes", c_to_n_ratio$phylum), "Bacteroidetes",as.character(c_to_n_ratio$class))
  c_to_n_ratio = c_to_n_ratio[,c(3,6,7,8)]
  
  # Filter by genera with significant assimilation (labeling) patterns
  keep_classes = c("Alphaproteobacteria","Bacteroidetes","Gammaproteobacteria")
  keep_genera = c("Amylibacter","Candidatus_Actinomarina","Flavivirga","Flavobacterium","Formosa","Gaetbulibacter","Gramella","Halocynthiibacter","Hellea","Kordia","Lacinutrix","Lentibacter","Lutibacter","Mesoflavibacter","OM182_clade","Pelagicola","Planktomarina","Polaribacter","Pseudooceanicola","Roseobacter","Roseovarius","SAR86_cluster","Sediminicola","Sulfitobacter","Tenacibaculum","Thalassobius","Unassigned_Alphaproteobacteria","Unassigned_Gammaproteobacteria","Unassigned_Methylophilales","Unassigned_Rhodobacterales","Unassigned_Actinomycetales","Unassigned_Alphaproteobacteria","Unassigned_Alteromonadales","Unassigned_Bacteroidetes","Unassigned_Candidatus_Actinomarinales","Unassigned_Cellvibrionales","Unassigned_Crocinitomicaceae","Unassigned_Cryomorphaceae","Unassigned_Euryarchaeota","Unassigned_Flavobacteriaceae","Unassigned_Flavobacteriales","Unassigned_Porticoccaceae","Unassigned_Pseudohongiella","Unassigned_Rhodobacteraceae","Unassigned_Rhodobacterales","Unassigned_Rhodospirillaceae","Winogradskyella")
  c_to_n_ratio = c_to_n_ratio[c_to_n_ratio$genus %in% keep_genera,]
  c_to_n_ratio = c_to_n_ratio[c_to_n_ratio$class %in% keep_classes,]
  
  # Aggregate into genera
  c_to_n_ratio_agg = aggregate(. ~ class + genus, c_to_n_ratio, mean)
  
  fig4c_plot = ggplot(c_to_n_ratio_agg, aes(x=aa_c_to_n,y=gc,fill=class)) +
    geom_point(shape=21,size=2.5,color="black",stroke=0.5) +
    scale_fill_manual(values = c("Alphaproteobacteria"='#2171B5', "Bacteroidetes" ='gray50', "Gammaproteobacteria" ='darkorange')) +
    theme_minimal() +
    xlab("Protein C:N") +
    ylab("DNA GC%") +
    theme(panel.border = element_rect(colour = "black", fill=NA),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.title = element_text(size=12,color="black"),
          axis.text.x=element_text(size=12,color="black"), 
          axis.text.y=element_text(size=12,color="black"))
  
  ggsave("fig4c.svg", fig4c_plot, device = "svg",scale = 1, width = 4.5, height = 2.5, units = "in", dpi = 600)
  
  
  
}

#####################################
###### SUPPLEMENTAL FIGURE 1 ######## 
#####################################
if(plot_sfig1 == "T"){
  
  # Read in data
  all_data = read.csv("all_protein_data.csv", header=T, sep=",", stringsAsFactors = F, comment.char = "")
  
  # Aggregate rows
  sample_by_enrichment_df = aggregate(. ~ species + sample, data=all_data[,c(1,11,12:ncol(all_data))], FUN=sum)
  #sample_by_enrichment_df = sample_by_enrichment_df[!grepl("^t",sample_by_enrichment_df$sample),]
  
  # Aggregate columns
  sample_by_enrichment_df$sum = rowSums(sample_by_enrichment_df[,c(3:4)])
  unlabeled_proteins = sample_by_enrichment_df[,c(1:2,ncol(sample_by_enrichment_df))]
  
  # Unmelt data to create relative abundance matrix
  unlabeled_proteins_unmelted = dcast(data = unlabeled_proteins,formula = sample~species,fun.aggregate = sum,value.var = "sum")
  rownames(unlabeled_proteins_unmelted) = unlabeled_proteins_unmelted$sample
  unlabeled_proteins_unmelted$sample = NULL
  unlabeled_proteins_unmelted = as.data.frame(t(unlabeled_proteins_unmelted))
  cols = as.character(colnames(unlabeled_proteins_unmelted))
  rows = as.character(rownames(unlabeled_proteins_unmelted))
  unlabeled_proteins_unmelted = unlabeled_proteins_unmelted %>% mutate_all(as.character)
  unlabeled_proteins_unmelted = sapply(unlabeled_proteins_unmelted, as.numeric)
  unlabeled_proteins_unmelted = sweep(unlabeled_proteins_unmelted, 2, colSums(unlabeled_proteins_unmelted), '/')
  rownames(unlabeled_proteins_unmelted) = rows
  colnames(unlabeled_proteins_unmelted) = cols
  
  # Calculate distance matrix
  samples_dist = vegdist(t(unlabeled_proteins_unmelted), method = "bray")
  
  # Calculate PCoA on distance matrix
  pcoa = pcoa(samples_dist)
  axis_1_variance = round(pcoa$values$Relative_eig[1]*100,digits = 2)
  axis_2_variance = round(pcoa$values$Relative_eig[2]*100, digits = 2)
  pcoa_vectors = data.frame(pcoa$vectors[,1:2])
  vectors_to_plot = pcoa_vectors
  colnames(vectors_to_plot) = c("Axis.1","Axis.2")
  vectors_to_plot$sample = rownames(vectors_to_plot)
  vectors_to_plot$treatment = substr(vectors_to_plot$sample,1,nchar(as.vector(vectors_to_plot$sample))-1)
  
  sfig1_plot = ggplot(data=vectors_to_plot, aes(x=Axis.1, y = Axis.2,fill=treatment, color=treatment)) +
    geom_hline(yintercept = 0, lty="dotted")+
    geom_vline(xintercept = 0, lty="dotted")+
    geom_polygon(size=0.5,alpha=0.2) +
    geom_point(shape=21, size= 2.5) +
    scale_color_manual(name="",guide = guide_legend(reverse = TRUE),values = c("CEx"="gray80","CLy"="gray60","DEx"="gray40","DLy"="gray20","t0"="gray90","t1"="black")) +
    scale_fill_manual(name="",guide = guide_legend(reverse = TRUE),values = c("CEx"="gray80","CLy"="gray60","DEx"="gray40","DLy"="gray20","t0"="gray90","t1"="black")) +
    xlab(paste("PCoA 1 (",axis_1_variance,"%)",sep="")) + 
    ylab(paste("PCoA 2 (",axis_2_variance,"%)",sep="")) + 
    xlim(c(-0.3,0.4)) +
    ylim(c(-0.2,0.14)) +
    theme(panel.border = element_rect(colour = "black", fill=NA),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.title = element_text(size=12,color="black"),
          axis.text.x=element_text(size=12,color="black"), 
          axis.text.y=element_text(size=12,color="black"))
  
  ggsave("sfig1.svg", sfig1_plot, device = "svg",scale = 1, width = 3.5, height = 2.5, units = "in", dpi = 600)
  
}

#####################################
###### SUPPLEMENTAL FIGURE 8 ######## 
#####################################
if(plot_sfig2 == "T"){
  
  # Read in data and filter
  all_data = read.csv("all_protein_data.csv", header=T, sep=",", stringsAsFactors = F, comment.char = "")
  all_data$class = ifelse(grepl("Bacteroidetes", all_data$phylum), "Bacteroidetes",as.character(all_data$class))
  keep_classes = c("Alphaproteobacteria","Bacteroidetes","Gammaproteobacteria")
  keep_genera = c("Amylibacter","Candidatus_Actinomarina","Flavivirga","Flavobacterium","Formosa","Gaetbulibacter","Gramella","Halocynthiibacter","Hellea","Kordia","Lacinutrix","Lentibacter","Lutibacter","Mesoflavibacter","OM182_clade","Pelagicola","Planktomarina","Polaribacter","Pseudooceanicola","Roseobacter","Roseovarius","SAR86_cluster","Sediminicola","Sphingobacterium","Sulfitobacter","Tenacibaculum","Thalassobius","Unassigned_Alphaproteobacteria","Unassigned_Gammaproteobacteria","Unassigned_Methylophilales","Unassigned_Rhodobacterales","Unassigned_Actinomycetales","Unassigned_Alphaproteobacteria","Unassigned_Alteromonadales","Unassigned_Bacteroidetes","Unassigned_Candidatus_Actinomarinales","Unassigned_Cellvibrionales","Unassigned_Crocinitomicaceae","Unassigned_Cryomorphaceae","Unassigned_Euryarchaeota","Unassigned_Flavobacteriaceae","Unassigned_Flavobacteriales","Unassigned_Porticoccaceae","Unassigned_Pseudohongiella","Unassigned_Rhodobacteraceae","Unassigned_Rhodobacterales","Unassigned_Rhodospirillaceae","Winogradskyella")
  all_data = all_data[all_data$class %in% keep_classes,]
  all_data = all_data[all_data$genus %in% keep_genera,]
  
  # Remove extraneous information
  sample_by_class_by_enrichment_df = aggregate(. ~ sample + class + genus, data=all_data[,c(1,7,10,12:ncol(all_data))], FUN=sum)
  sample_by_class_by_enrichment_df = sample_by_class_by_enrichment_df[!grepl("^t",sample_by_class_by_enrichment_df$sample),]
  
  # Define labeled (>= 2%) and unlabeled (< 2%) fractions 
  samples = sample_by_class_by_enrichment_df$sample
  classes = sample_by_class_by_enrichment_df$class
  genera = sample_by_class_by_enrichment_df$genus
  all_proteins = sample_by_class_by_enrichment_df[,4:ncol(sample_by_class_by_enrichment_df)]
  unlabeled_proteins = sample_by_class_by_enrichment_df[,4:5]
  labeled_proteins = sample_by_class_by_enrichment_df[,6:ncol(sample_by_class_by_enrichment_df)]
  
  # Calculate ENRV
  enrv_multiplier = as.numeric(2:100)
  enrv = data.frame(mapply(`*`,labeled_proteins,enrv_multiplier))
  enrv = rowSums(enrv)/rowSums(labeled_proteins)
  
  # Calculate LF
  lf = (rowSums(labeled_proteins)/rowSums(all_proteins))*100
  
  # Create enrv summary
  enrv_summary = data.frame(samples, paste(classes,genera,sep="__"),enrv)
  colnames(enrv_summary) = c("sample","taxon","enrv")
  
  # Unmelt to normalize and fill in NAs
  enrv_summary_unmelted = dcast(enrv_summary, taxon ~ sample, value.var = "enrv",sum)
  rownames(enrv_summary_unmelted) = enrv_summary_unmelted$taxon
  enrv_summary_unmelted$taxon = NULL
  enrv_summary_unmelted[is.na(enrv_summary_unmelted)] = 0
  enrv_summary_unmelted = sweep(enrv_summary_unmelted,2,colSums(enrv_summary_unmelted),`/`)
  
  # Create distance matrix based on enrv values
  enrv_summary_dist = vegdist(enrv_summary_unmelted, method = "bray")
  
  # Create taxon dendrogram based on distance matrix
  enrv_summary_hclust = hclust(enrv_summary_dist)
  enrv_summary_hclust_labs = as.data.frame(enrv_summary_hclust$labels)
  enrv_summary_hclust_labs$id <- as.integer(row.names(enrv_summary_hclust_labs)) 
  enrv_summary_hclust_labs = enrv_summary_hclust_labs[order(match(enrv_summary_hclust_labs$id, enrv_summary_hclust$order)), , drop = FALSE]
  taxaorder = as.character(enrv_summary_hclust_labs[,1])
  
  # Calculate Z-score
  enrv_summary = as.data.frame(melt(t(as.matrix(enrv_summary_unmelted))))
  colnames(enrv_summary) = c("sample","taxon","enrv")
  enrv_summary$zscore = ave(enrv_summary$enrv, enrv_summary$taxon, FUN = scale)
  enrv_summary$sample = substr(enrv_summary$sample,1,nchar(as.vector(enrv_summary$sample))-1)
  enrv_summary_agg = aggregate(zscore ~ sample + taxon, enrv_summary, mean)
  
  # Order rows and columns for heatmap
  enrv_summary_agg$taxon = reorder.factor(enrv_summary_agg$taxon, new.order=taxaorder)
  enrv_summary_agg$sample = factor(enrv_summary_agg$sample, levels = c("DLy","DEx","CLy","CEx"))
  
  # Plot
  sfig2_plot = ggplot(data=enrv_summary_agg,aes(x=sample,y=taxon,fill=zscore)) +
    geom_tile(color="black", size=0.1) +
    xlab("") +
    ylab("") +
    #scale_fill_viridis(name="Norm Spec Cts")+
    scale_fill_gradient(low = "navyblue",high = "yellow",name="Per-taxon\nZ-score of ENRV")+
    theme_minimal()
  
  ggsave("sfig2.svg", sfig2_plot, device = "svg",scale = 1, width = 6, height = 8, units = "in", dpi = 600)
  
}

#####################################
###### SUPPLEMENTAL FIGURE 4 ######## 
#####################################
if(plot_sfig4 == "T"){
  
  # Read in data
  all_data = read.csv("all_protein_data.csv", header=T, sep=",", stringsAsFactors = F, comment.char = "")
  all_data$class = ifelse(grepl("Bacteroidetes", all_data$phylum), "Bacteroidetes",as.character(all_data$class))
  
  # Aggregate rows
  sample_by_enrichment_df = aggregate(. ~ class + order + sample, data=all_data[,c(1,7,8,12:ncol(all_data))], FUN=sum)
  sample_by_enrichment_df = sample_by_enrichment_df[!grepl("^t",sample_by_enrichment_df$sample),]
  
  # Aggregate columns
  sample_by_enrichment_df$sum = rowSums(sample_by_enrichment_df[,4:5])
  unlabeled_proteins = sample_by_enrichment_df[,c(1:3,ncol(sample_by_enrichment_df))]
  unlabeled_proteins$taxon = paste(unlabeled_proteins$class,unlabeled_proteins$order,sep="__")
  unlabeled_proteins$sample = substr(unlabeled_proteins$sample,1,nchar(as.vector(unlabeled_proteins$sample))-1)
  unlabeled_proteins = aggregate(sum ~ taxon + sample, unlabeled_proteins, FUN=mean)
  
  # Unmelt data to create relative abundance matrix
  unlabeled_proteins_unmelted = dcast(data = unlabeled_proteins,formula = sample~taxon,fun.aggregate = sum,value.var = "sum")
  rownames(unlabeled_proteins_unmelted) = unlabeled_proteins_unmelted$sample
  unlabeled_proteins_unmelted$sample = NULL
  unlabeled_proteins_unmelted = as.data.frame(t(unlabeled_proteins_unmelted))
  cols = as.character(colnames(unlabeled_proteins_unmelted))
  rows = as.character(rownames(unlabeled_proteins_unmelted))
  unlabeled_proteins_unmelted = unlabeled_proteins_unmelted %>% mutate_all(as.character)
  unlabeled_proteins_unmelted = sapply(unlabeled_proteins_unmelted, as.numeric)
  unlabeled_proteins_unmelted = sweep(unlabeled_proteins_unmelted, 2, colSums(unlabeled_proteins_unmelted), '/')
  rownames(unlabeled_proteins_unmelted) = rows
  colnames(unlabeled_proteins_unmelted) = cols
  unlabeled_proteins_unmelted = as.data.frame(unlabeled_proteins_unmelted)
  unlabeled_proteins_unmelted$taxon = rownames(unlabeled_proteins_unmelted)
  unlabeled_proteins_unmelted = cSplit(unlabeled_proteins_unmelted, "taxon", "__")
  
  keep_taxa = c("Rhodobacterales","Pelagibacterales","Rhodospirillales","SAR116_cluster","Rhizobiales","Unassigned_Alphaproteobacteria", ## Alphas
                "Flavobacteriales","Sphingobacteriales","Bacteroidales","Cytophagales", "Unassigned_Bacteroidetes",  ## Bacts
                "SAR86_cluster", "Alteromonadales","Cellvibrionales","Oceanospirillales", "Unassigned_Gammaproteobacteria", ## Gammas
                "Unassigned_Euryarchaeota",  ## Archaea
                "Actinobacteria",  ## Actinos
                "Burkholderiales","Methylophilales",  ## Betas
                "Unclassified",   ## Unclassified
                "Other")
  
  # Add "Other" category
  unlabeled_proteins_unmelted_subset = data.frame(unlabeled_proteins_unmelted[unlabeled_proteins_unmelted$taxon_2 %in% keep_taxa,])
  unlabeled_proteins_unmelted_subset = transform(unlabeled_proteins_unmelted_subset, taxon_1 = as.character(taxon_1),taxon_2 = as.character(taxon_2),CEx = as.numeric(CEx),CLy = as.numeric(CLy),DEx = as.numeric(DEx),DLy = as.numeric(DLy))
  other = data.frame(rbind(1-colSums(unlabeled_proteins_unmelted_subset[,1:4])),"Other","Other")
  colnames(other)[5:6] = c("taxon_1","taxon_2")
  unlabeled_proteins_unmelted_subset = rbind(unlabeled_proteins_unmelted_subset,other)
  unlabeled_proteins_unmelted_subset$taxon_1 = ifelse(!grepl("Alphaproteobacteria|Gammaproteobacteria|Bacteroidetes", unlabeled_proteins_unmelted_subset$taxon_1),"Other",unlabeled_proteins_unmelted_subset$taxon_1)
  unlabeled_proteins_unmelted_subset$sum = rowSums(unlabeled_proteins_unmelted_subset[,1:4])
  unlabeled_proteins_unmelted_subset = unlabeled_proteins_unmelted_subset[order(-unlabeled_proteins_unmelted_subset$sum),]
  unlabeled_proteins_unmelted_subset$sum = NULL
  
  # Create barplot aes
  barplot_data = split(unlabeled_proteins_unmelted_subset,f = unlabeled_proteins_unmelted_subset$taxon_1)
  barplot_data$Alphaproteobacteria$colors = brewer.pal(length(as.character(unique(barplot_data$Alphaproteobacteria$taxon_2))), "Blues")
  barplot_data$Bacteroidetes$colors = brewer.pal(length(as.character(unique(barplot_data$Bacteroidetes$taxon_2))), "Greys")
  barplot_data$Gammaproteobacteria$colors = brewer.pal(length(as.character(unique(barplot_data$Gammaproteobacteria$taxon_2))), "Oranges")
  barplot_data$Other$colors = brewer.pal(length(as.character(unique(barplot_data$Other$taxon_2))), "Purples")
  
  # Melt data into useable form
  barplot_data$Alphaproteobacteria$taxon_2 = factor(barplot_data$Alphaproteobacteria$taxon_2, levels = barplot_data$Alphaproteobacteria$taxon_2)
  barplot_data$Alphaproteobacteria = melt(barplot_data$Alphaproteobacteria)
  barplot_data$Alphaproteobacteria$variable = factor(barplot_data$Alphaproteobacteria$variable, levels = c("DLy","DEx","CLy","CEx"))
  barplot_data$Bacteroidetes$taxon_2 = factor(barplot_data$Bacteroidetes$taxon_2, levels = barplot_data$Bacteroidetes$taxon_2)
  barplot_data$Bacteroidetes = melt(barplot_data$Bacteroidetes)
  barplot_data$Bacteroidetes$variable = factor(barplot_data$Bacteroidetes$variable, levels = c("DLy","DEx","CLy","CEx"))
  barplot_data$Gammaproteobacteria$taxon_2 = factor(barplot_data$Gammaproteobacteria$taxon_2, levels = barplot_data$Gammaproteobacteria$taxon_2)
  barplot_data$Gammaproteobacteria = melt(barplot_data$Gammaproteobacteria)
  barplot_data$Gammaproteobacteria$variable = factor(barplot_data$Gammaproteobacteria$variable, levels = c("DLy","DEx","CLy","CEx"))
  barplot_data$Other$taxon_2 = factor(barplot_data$Other$taxon_2, levels = barplot_data$Other$taxon_2)
  barplot_data$Other = melt(barplot_data$Other)
  barplot_data$Other$variable = factor(barplot_data$Other$variable, levels = c("DLy","DEx","CLy","CEx"))
  
  
  barplot1 = ggplot(data=barplot_data$Alphaproteobacteria, aes(x=variable, y=value, group=rev(taxon_2), fill=taxon_2)) +
    geom_bar(show.legend = T,stat="identity", position = position_stack(reverse=TRUE),colour = "gray40", width=0.9) +
    scale_y_continuous(expand = c(0.05,0), limits=c(0,0.6)) +
    scale_fill_manual(guide = "legend", values = barplot_data$Alphaproteobacteria$colors) +
    theme_minimal() +
    xlab("") +
    ylab("") +
    guides(fill=guide_legend(title="New Legend Title")) +
    theme(axis.text.y = element_text(size=10, colour="black"),
          axis.text.x = element_text(size = 10, colour = "black"), 
          axis.title.x=element_blank(),
          legend.text = element_text(size=10),
          legend.key.size = unit(0.6,"cm"),
          legend.title = element_blank())
  
  # Define other barplots based on the original
  barplot2 = barplot1 %+% barplot_data$Bacteroidetes + scale_fill_manual(guide = "legend", values = barplot_data$Bacteroidetes$colors)
  barplot3 = barplot1 %+% barplot_data$Gammaproteobacteria + scale_fill_manual(guide = "legend", values = barplot_data$Gammaproteobacteria$colors)
  barplot4 = barplot1 %+% barplot_data$Other + scale_fill_manual(guide = "legend", values = barplot_data$Other$colors)
  
  # Extract legends
  barplot1_legend = get_legend(barplot1)
  barplot2_legend = get_legend(barplot2)
  barplot3_legend = get_legend(barplot3)
  barplot4_legend = get_legend(barplot4)
  
  # Remove legends for first facet
  barplot1 = barplot1 + theme(legend.position = "none") 
  barplot2 = barplot2 + theme(legend.position = "none") 
  barplot3 = barplot3 + theme(legend.position = "none") 
  barplot4 = barplot4 + theme(legend.position = "none") 
  
  # Create barplot facet and legend facet
  all_barplots = grid.arrange(barplot1, barplot2, barplot3, barplot4,ncol = 1)
  all_barplot_legends = grid.arrange(barplot1_legend, barplot2_legend, barplot3_legend, barplot4_legend,ncol = 1)
  
  # Stick together bars and legends
  sfig4_plot = grid.arrange(all_barplots,all_barplot_legends,ncol=2)
  
  # Print
  ggsave("sfig4.svg", sfig4_plot, device = "svg",scale = 1, width = 5, height = 6, units = "in", dpi = 600)
  
}

#####################################
###### SUPPLEMENTAL FIGURE 6 ######## 
#####################################
if(plot_sfig6 == "T"){
  
  # Read in data
  all_data = read.csv("all_protein_data.csv", header=T, sep=",", stringsAsFactors = F, comment.char = "")
  all_data$class = ifelse(grepl("Bacteroidetes", all_data$phylum), "Bacteroidetes",as.character(all_data$class))

  # Aggregate rows
  sample_by_enrichment_df = aggregate(. ~  order + sample, data=all_data[,c(1,8,12:ncol(all_data))], FUN=sum)
  keep_orders = c("Rhodobacterales","Pelagibacterales","Flavobacteriales","Rhodospirillales","Cellvibrionales","SAR86_cluster","Oceanospirillales", "Nitrosomonadales","Marine_Group_II_Euryarchaeota", "Acidimicrobiales")
  sample_by_enrichment_df = sample_by_enrichment_df[sample_by_enrichment_df$order %in% keep_orders,]

  # Create unlabeled x labeled square matrix
  unlabeled_proteins = rowSums(sample_by_enrichment_df[,3:4])
  labeled_proteins = rowSums(sample_by_enrichment_df[,5:ncol(sample_by_enrichment_df)])
  deviation_calculations = cbind(sample_by_enrichment_df[,1:2],unlabeled_proteins,labeled_proteins)
  deviation_calculations$total_proteins = deviation_calculations$unlabeled_proteins + deviation_calculations$labeled_proteins
  deviation_calculations$sample = substr(deviation_calculations$sample,1,nchar(as.vector(deviation_calculations$sample))-1) 
  deviation_calculations_agg = aggregate(. ~ sample + order, deviation_calculations,FUN=mean)

  # Calculate relative abundances of groups
  deviation_calculations_relative = deviation_calculations_agg %>% group_by(sample) %>% mutate(unlabeled_proteins_relative = (unlabeled_proteins / sum(unlabeled_proteins)*100))
  deviation_calculations_relative = deviation_calculations_relative %>% group_by(sample) %>% mutate(labeled_proteins_relative = (labeled_proteins / sum(labeled_proteins)*100))
  deviation_calculations_relative = deviation_calculations_relative %>% group_by(sample) %>% mutate(total_proteins_relative = (total_proteins / sum(total_proteins)*100))
  
  # Get subsets to plot
  total_t1_data = deviation_calculations_relative[deviation_calculations_relative$sample == "t1",][,c(1,2,8)] 
  total_treatment_data = deviation_calculations_relative[!grepl("^t",deviation_calculations_relative$sample),][,c(1,2,7,8)] 
  
  # Merge into a single plotting df
  plotting_data = data.frame(merge(total_t1_data, total_treatment_data, by="order",all.x = T))[,c(1,4,3,5,6)]
  colnames(plotting_data) = c("taxon","sample","control_relabund","labeled_treatment_relabund","total_treatment_relabund") 
  plotting_data$taxon = factor(plotting_data$taxon, levels = keep_orders)
  plotting_data$sample = factor(plotting_data$sample, levels = c("DLy","DEx","CLy","CEx"))
  
  scatterplot1 = ggplot(plotting_data, aes(x=log10(control_relabund+1), y=log10(total_treatment_relabund+1), color=sample, shape=taxon)) +
    geom_point(size=2, show.legend = T, stroke=1) +
    scale_shape_manual(values=rev(1:nlevels(plotting_data$taxon))) +
    scale_color_manual(values = c("gray20","gray40","gray60","gray90")) +
    scale_y_continuous(limits = c(0,2)) +
    scale_x_continuous(limits = c(0,2)) +
    geom_abline(slope = 1, lty="dashed", size=0.25) +
    ylab("% of total substrate proteomes (log)") +
    xlab("% of T1 proteome (log)") +
    theme_classic() +
    theme(panel.spacing = unit(2, "lines"),
          aspect.ratio = 1,
          axis.line = element_line(),
          axis.text.x = element_text(color="black"),
          axis.text.y = element_text(color="black"),
          strip.text = element_blank(),
          strip.background = element_blank())
  
  scatterplot2 = ggplot(plotting_data, aes(x=log10(total_treatment_relabund+1), y=log10(labeled_treatment_relabund+1), color=sample, shape=taxon)) +
    geom_point(size=2, show.legend = T, stroke=1) +
    scale_shape_manual(values=rev(1:nlevels(plotting_data$taxon))) +
    scale_color_manual(values = c("gray20","gray40","gray60","gray90")) +
    scale_y_continuous(limits = c(0,2)) +
    scale_x_continuous(limits = c(0,2)) +
    geom_abline(slope = 1, lty="dashed", size=0.25) +
    ylab("% of labeled substrate proteomes (log)") +
    xlab("% of total substrate proteomes (log)") +
    theme_classic() +
    theme(panel.spacing = unit(2, "lines"),
          aspect.ratio = 1,
          axis.line = element_line(),
          axis.text.x = element_text(color="black"),
          axis.text.y = element_text(color="black"),
          strip.text = element_blank(),
          strip.background = element_blank())
  
  # Extract legend
  scatterplot_legend = get_legend(scatterplot1)
  
  # Remove legends for first facet
  scatterplot1 = scatterplot1 + theme(legend.position = "none")
  scatterplot2 = scatterplot2 + theme(legend.position = "none")
  
  stacked_scatterplots = grid.arrange(scatterplot1, scatterplot2,ncol=1)
  sfig6_plot = grid.arrange(stacked_scatterplots,scatterplot_legend, ncol=2)
  
  ggsave("sfig6.svg", sfig6_plot, device = "svg",scale = 1, width = 5, height = 6, units = "in", dpi = 600)
  
  
}

#####################################
###### SUPPLEMENTAL FIGURE 7 ######## 
#####################################
if(plot_sfig7 == "T"){

  # Read in data
  all_data = read.csv("all_protein_data.csv", header=T, sep=",", stringsAsFactors = F, comment.char = "")
  all_data$class = ifelse(grepl("Bacteroidetes", all_data$phylum), "Bacteroidetes",as.character(all_data$class))
  
  # Define plotting variables based on abundant enzymes/pathways/functions
  all_data$cog_desc = ifelse(grepl("^Ribosomal_protein",all_data$cog_desc),"Ribosomal_proteins",as.character(all_data$cog_desc))
  all_data$cog_desc = ifelse(grepl("polymerase",all_data$cog_desc),"DNA/RNA_polymerases",as.character(all_data$cog_desc))
  all_data$cog_desc = ifelse(grepl("^F0F1",all_data$cog_desc),"ATPases",as.character(all_data$cog_desc))
  all_data$cog_desc = ifelse(grepl("vacuolar-type",all_data$cog_desc),"ATPases",as.character(all_data$cog_desc))
  all_data$cog_desc = ifelse(grepl("^MoxR",all_data$cog_desc),"ATPases",as.character(all_data$cog_desc))
  all_data$cog_desc = ifelse(grepl("ranscription",all_data$cog_desc),"Transcription/Translation_factors",as.character(all_data$cog_desc))
  all_data$cog_desc = ifelse(grepl("ranslation",all_data$cog_desc),"Transcription/Translation_factors",as.character(all_data$cog_desc))
  all_data$cog_desc = ifelse(grepl("Polyribonucleotide",all_data$cog_desc),"Transcription/Translation_factors",as.character(all_data$cog_desc))
  all_data$cog_desc = ifelse(grepl("Biopolymer_transport_protein",all_data$cog_desc),"Biopolymer_transport",as.character(all_data$cog_desc))
  all_data$cog_desc = ifelse(grepl("Flagellar_motor_protein",all_data$cog_desc),"Biopolymer_transport",as.character(all_data$cog_desc))
  all_data$cog_desc = ifelse(grepl("Periplasmic_protein_involved_in_polysaccharide_export",all_data$cog_desc),"Biopolymer_transport",as.character(all_data$cog_desc))
  all_data$cog_desc = ifelse(grepl("Poly(3-hydroxybutyrate)_depolymerase",all_data$cog_desc),"Biopolymer_transport",as.character(all_data$cog_desc))
  all_data$cog_desc = ifelse(grepl("permidine",all_data$cog_desc),"Polyamine_transport",as.character(all_data$cog_desc))
  all_data$cog_desc = ifelse(grepl("ABC-type",all_data$cog_desc),"ABC_transport",as.character(all_data$cog_desc))
  all_data$cog_desc = ifelse(grepl("ABC_trans",all_data$cog_desc),"ABC_transport",as.character(all_data$cog_desc))
  all_data$cog_desc = ifelse(grepl("TRAP-type",all_data$cog_desc),"TRAP_transport",as.character(all_data$cog_desc))
  all_data$cog_desc = ifelse(grepl("Outer_membrane_receptor_proteins_",all_data$cog_desc),"TonB_transport",as.character(all_data$cog_desc))
  all_data$cog_desc = ifelse(grepl("TonB",all_data$cog_desc),"TonB_transport",as.character(all_data$cog_desc))
  all_data$cog_desc = ifelse(grepl("Outer_membrane_cobalamin_receptor_protein",all_data$cog_desc),"TonB_transport",as.character(all_data$cog_desc))
  all_data$cog_desc = ifelse(grepl("Outer_membrane_receptor_for_ferrienterochelin_and_colicins",all_data$cog_desc),"TonB_transport",as.character(all_data$cog_desc))
  all_data$cog_desc = ifelse(grepl("Riboso",all_data$cog_desc),"Ribosomal_proteins",as.character(all_data$cog_desc))
  all_data$cog_desc = gsub('_system_periplasmic_component', '', as.character(all_data$cog_desc))
  all_data$cog_desc = gsub('_systems_periplasmic_component', '', as.character(all_data$cog_desc))
  all_data$cog_desc = gsub('/domain', '', as.character(all_data$cog_desc))
  
  # Aggregate rows
  sample_by_enrichment_df = aggregate(. ~ cog_desc + class + sample, data=all_data[,c(1,4,7,12:ncol(all_data))], FUN=sum)
  keep_classes = c("Alphaproteobacteria","Bacteroidetes","Gammaproteobacteria")
  keep_functions = c("Other","ABC_transport","TRAP_transport","TonB_transport","Polyamine_transport","Biopolymer_transport","DNA/RNA_polymerases", "Nucleoside_diphosphate_kinase", "Ribosomal_proteins", "Transcription/Translation_factors")
  sample_by_enrichment_df = sample_by_enrichment_df[sample_by_enrichment_df$class %in% keep_classes,]
  sample_by_enrichment_df$cog_desc = ifelse(!(sample_by_enrichment_df$cog_desc %in% keep_functions),"Other",sample_by_enrichment_df$cog_desc)
  sample_by_enrichment_df = sample_by_enrichment_df[sample_by_enrichment_df$cog_desc %in% keep_functions,]
  sample_by_enrichment_df = aggregate(. ~ cog_desc + class + sample,sample_by_enrichment_df,FUN=sum)
  
  # Create unlabeled x labeled square matrix
  unlabeled_proteins = rowSums(sample_by_enrichment_df[,4:5])
  labeled_proteins = rowSums(sample_by_enrichment_df[,6:ncol(sample_by_enrichment_df)])
  deviation_calculations = cbind(sample_by_enrichment_df[,1:3],unlabeled_proteins,labeled_proteins)
  deviation_calculations$total_proteins = deviation_calculations$unlabeled_proteins + deviation_calculations$labeled_proteins
  deviation_calculations_agg = aggregate(. ~ sample + class + cog_desc, deviation_calculations,FUN=mean)
  
  # Calculate relative abundances of groups
  deviation_calculations_relative = deviation_calculations_agg %>% group_by(sample, class) %>% mutate(unlabeled_proteins_relative = (unlabeled_proteins / sum(unlabeled_proteins)*100))
  deviation_calculations_relative = deviation_calculations_relative %>% group_by(sample, class) %>% mutate(labeled_proteins_relative = (labeled_proteins / sum(labeled_proteins)*100))
  deviation_calculations_relative = deviation_calculations_relative %>% group_by(sample, class) %>% mutate(total_proteins_relative = (total_proteins / sum(total_proteins)*100))
  
  # Define plotted data
  deviation_calculations_relative$relative_data_to_plot = ifelse(grepl("^t",deviation_calculations_relative$sample),deviation_calculations_relative$total_proteins_relative,deviation_calculations_relative$labeled_proteins_relative)
  deviation_calculations_relative$library = substr(deviation_calculations_relative$sample,1,nchar(as.vector(deviation_calculations_relative$sample))-1)
  deviation_calculations_relative$library = factor(deviation_calculations_relative$library, levels = c("t0","t1","DLy","DEx","CLy","CEx"))
  deviation_calculations_relative$cog_desc = factor(deviation_calculations_relative$cog_desc, levels = keep_functions)
  
  sfig7_plot = ggplot(deviation_calculations_relative, aes(x=cog_desc, y=relative_data_to_plot,color=library,fill =class)) +
    geom_violin(scale = "width",draw_quantiles = 0.5) +
    scale_color_manual(values = c("black","gray90","gray20","gray40","gray60","gray80")) +
    scale_fill_manual(values = c("Alphaproteobacteria"='#2171B5', "Bacteroidetes" ='gray50', "Gammaproteobacteria" ='darkorange')) +
    facet_wrap(~class, scales = "fixed", ncol = 1) +
    theme_minimal() +
    xlab("") +
    ylab("Relative proteome abundance (%)") +
    theme(axis.text.x =  element_text(color = "black",angle=45, hjust=1,size=5),
          axis.text.y = element_text(size=10, colour="black"),
          panel.border = element_rect(colour = "black", fill=NA))
    
  ggsave("sfig7.svg", sfig7_plot, device = "svg",scale = 1, width = 6, height = 6, units = "in", dpi = 600)
  
}








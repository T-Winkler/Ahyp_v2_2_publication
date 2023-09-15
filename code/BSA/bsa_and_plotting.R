
setwd("/home/tom/Documents/projects/Ahyp_v2_2_publication/")

library(QTLseqr)
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())
library(gggenes)
library(ggtranscript)
library(ape)
library(patchwork)


bsa_analysis <- function(rawData,HighBulk,LowBulk,Chroms,nhigh,nlow){
  df <- importFromGATK(file = rawData,
                       highBulk = HighBulk,
                       lowBulk = LowBulk,
                       chromList = Chroms)
  
  df <- df %>% select(CHROM,
                      POS,
                      REF,
                      ALT,
                      AD_REF.LOW,
                      AD_ALT.LOW,
                      DP.LOW,
                      GQ.LOW,
                      PL.LOW,
                      SNPindex.LOW,
                      AD_REF.HIGH,
                      AD_ALT.HIGH,
                      DP.HIGH,
                      GQ.HIGH,
                      PL.HIGH,
                      SNPindex.HIGH,
                      REF_FRQ,
                      deltaSNP) %>%
    mutate(CHROM=as.factor(as.numeric(gsub("Scaffold_","",CHROM))),
           POS=as.numeric(POS)) %>%
    filter(REF!='*',ALT!='*')
  
  df_filt <-filterSNPs(SNPset = df,
                       refAlleleFreq = 0.2,
                       minTotalDepth = 50,
                       maxTotalDepth = 100,
                       minSampleDepth = 20,
                       minGQ = 99,
                       verbose = TRUE)
  
  df_filt <- runGprimeAnalysis(
    SNPset = df_filt,
    windowSize = 2e6,
    filter = 0.4,
    outlierFilter = "deltaSNP")
  df_filt <- runQTLseqAnalysis(
    SNPset = df_filt,
    windowSize = 2e6,
    popStruc = "RIL",
    bulkSize =  c(nhigh, nlow),
    replications = 10000,
    filter = 0.4,
    intervals = c(95, 99)
  )
  return(df_filt)
}

AM_00332_leaf_green_red <- bsa_analysis(rawData = 'data/BSA/wgs/vcf/bulk_snps05.table',
                                        HighBulk = "AM_00332_gl",
                                        LowBulk = "AM_00332_rl",
                                        Chroms = paste0(rep("Scaffold_",
                                                            16),1:16),
                                        nhigh=80,
                                        nlow=80)
AM_00331_flower_red_green <- bsa_analysis(rawData = 'data/BSA/wgs/vcf/bulk_snps05.table',
                                          HighBulk = "AM_00331_rf",
                                          LowBulk = "AM_00331_gf",
                                          Chroms = paste0(rep("Scaffold_",
                                                              16),1:16),
                                          nhigh = 68,
                                          nlow = 68)


# plot all results
# leaf
plotGresults <- function(Gresults,betalain_genes){
  qval <- Gresults %>% 
    filter(qvalue<=0.01) 
  #qval <- min(qval$Gprime)
  qval <- 3
  
  
  mG <- Gresults %>%
    filter(Gprime==max(Gresults$Gprime))
  
  p1 <- ggplot()+
    geom_line(data=Gresults,aes(POS/1e6,Gprime), size=2) +
    labs(x= 'Position (Mb)',y= "G' value") +
    scale_x_continuous(breaks = c(0,10,20,30))+
    geom_hline(data=data.frame(yint=qval),
               aes(yintercept =yint,
                   linetype ='dashed',
                   color=alpha('red',0.6)), 
               size=1.7)+
    facet_grid(.~CHROM, space = 'free_x',scales='free_x') +
    theme(panel.spacing.x=unit(0.25, "lines")) +
    ylim(0,10) +
    theme(strip.background = element_rect(fill = alpha('lightblue',0.2)),
          strip.text = element_text(size=30)) +
    theme(legend.position="none",
          axis.text.y = element_text(size=40),
          axis.title.y = element_text(size=40),
          #axis.title.x = element_blank(),
          axis.title.x = element_text(size=40),
          axis.text.x = element_text(size = 20),
          panel.spacing.x = unit(6, "mm"),
          axis.line = element_line(linewidth = 2),
          axis.ticks = element_line(linewidth = 1.5),
          axis.ticks.length = unit(.25, "cm")) +
    scale_x_continuous(guide = guide_axis(check.overlap = T))
    #geom_gene_arrow(data=betalain_genes, 
    #                aes(xmin = start/1e6, xmax = end/1e6, y = max(Gresults$Gprime), fill = type))
  
  #ggsave(outfile,p1,width = 18,height = 7,,bg='white')
  return(p1)
}
# flower
plotGresults1 <- function(Gresults,betalain_genes){
  qval <- Gresults %>% 
    filter(qvalue<=0.01) 
  #qval <- min(qval$Gprime)
  qval <- 3
  
  
  mG <- Gresults %>%
    filter(Gprime==max(Gresults$Gprime))
  
  p1 <- ggplot()+
    geom_line(data=Gresults,aes(POS/1e6,Gprime), size=2) +
    labs(x= 'Position (Mb)',y= "G' value") +
    #labs(x= '',y= "G' value") +
    scale_x_continuous(breaks = c(0,10,20,30)) +
    geom_hline(data=data.frame(yint=qval),
               aes(yintercept =yint,
                   linetype ='dashed',
                   color=alpha('red',0.6)), 
               size=1.7)+
    facet_grid(.~CHROM,space = 'free_x',scales='free_x') +
    theme(panel.spacing.x=unit(0.25, "lines")) +
    ylim(0,8) +
    theme( strip.background = element_rect(fill = alpha('lightblue',0.2)),
           strip.text = element_text(size=30)) +
    theme(legend.position="none",
          axis.text.y = element_text(size=40),
          axis.title = element_text(size=40),
          axis.text.x = element_text(size = 20),
          panel.spacing.x = unit(6, "mm"),
          axis.line = element_line(linewidth = 2),
          axis.ticks = element_line(linewidth = 1.5),
          axis.ticks.length = unit(.25, "cm")) +
    scale_x_continuous(guide = guide_axis(check.overlap = T))
    #geom_gene_arrow(data=betalain_genes, 
    #                aes(xmin = start/1e6, xmax = end/1e6, y = max(Gresults$Gprime), fill = type))
  
  #ggsave(outfile,p1,width = 18,height = 7,,bg='white')
  return(p1)
}

# plot individual chromosomes
# flower
plotGqtl <- function(Gresults,chr,genes){
  
  qval <- Gresults %>% 
    filter(qvalue<=0.01) 
  #qval <- min(qval$Gprime)
  qval <- 3
  my_qtl <- getQTLTable(SNPset = Gresults, alpha = 0.01,export = F)
  
  ggplot()+
    geom_line(data=filter(Gresults,CHROM==chr),aes(POS/1e6,Gprime),size=2) +
    labs(x= 'Position (Mb)',y= "G' value") +
    scale_x_continuous(breaks = c(0,10,20,30))+
    geom_hline(data=data.frame(yint=qval),
               aes(yintercept = yint, 
                   linetype ='dashed', 
                   color=alpha('red',0.6)),
               size = 2) +
    facet_grid(.~CHROM,space = 'free_x',scales='free_x') + 
    theme(panel.spacing.x=unit(0.25, "lines")) +
    ylim(0,10) +
    theme( strip.background = element_rect(fill = alpha('lightblue',0.2)),
           strip.text = element_text(size=30)) +
    theme(legend.position="none",
          axis.text.y = element_blank(),
          axis.text.x = element_text(size=20),
          #axis.title.x = element_blank(),
          axis.title.x = element_text(size=40),
          axis.title.y = element_text(color = "white", size = 35),
          #axis.title.y = element_blank(),
          axis.line = element_line(linewidth = 2),
          axis.ticks = element_line(linewidth = 1.5),
          axis.ticks.length = unit(.25, "cm")) +
    geom_gene_arrow(data=filter(genes, CHROM==chr, type == "transcript") %>% droplevels(), 
                    aes(xmin = start/1e6, xmax = end/1e6, y = 9.5, color = attributes), size=6) +
    scale_color_manual(values = c(alpha('red',0.6), "black", "black","grey","grey"))
}
# leaf
plotGqtl1 <- function(Gresults,chr,genes){
  
  qval <- Gresults %>% 
    filter(qvalue<=0.01) 
  #qval <- min(qval$Gprime)
  qval <- 3
  my_qtl <- getQTLTable(SNPset = Gresults, alpha = 0.01,export = F)
  
  p1 <- ggplot() +
    geom_line(data=filter(Gresults,CHROM==chr),aes(POS/1e6,Gprime),size=2) +
    labs(x= 'Position (Mb)',y= "G' value") +
    scale_x_continuous(breaks = c(0,5,10,15,20,30))+
    scale_y_continuous(breaks = c(0,2,4,8)) +
    geom_hline(data=data.frame(yint=qval),
               aes(yintercept =yint, linetype ='dashed', color=alpha('red',0.6)),
               size = 2) +
    geom_vline(aes(xintercept = 5231549/1e6), color = "black") +
    geom_vline(aes(xintercept = 5305973/1e6), color = "black") +
    facet_grid(.~CHROM,space = 'free_x',scales='free_x') +
    ylim(0,8) +
    theme(panel.spacing.x=unit(0.25, "lines")) +
    # theme( strip.background = element_rect(fill = alpha('lightblue',0.2)),
    #        strip.text = element_text(size=30)) +
    theme( strip.background = element_blank(),
           strip.text = element_blank()) +
    theme(legend.position="none",
          axis.text.y = element_text(size = 30),
          #plot.margin = unit(c(0, 0, 0, 0), "cm"),
          axis.text.x = element_text(size=20),
          axis.title.x = element_text(size=40),
          #axis.title.y = element_blank(),
          axis.title.y = element_text(color = "black", size = 35),
          axis.line = element_line(linewidth = 2),
          axis.ticks = element_line(linewidth = 1.5),
          axis.ticks.length = unit(.25, "cm"))
  
  
  p2 <- ggplot() +
    geom_gene_arrow(data=filter(genes, 
                                CHROM==chr, 
                                type == "transcript",
                                attributes == "ID=AHp023147.1;geneID=AHp023147" | attributes == "ID=AHp023148.1;geneID=AHp023148") %>% droplevels(),
                    aes(xmin = start, 
                        xmax = end, 
                        y = "chr16", 
                        fill = attributes, 
                        forward = c(F,T)),
                    size = 1.5,
                    color = "black",
                    arrowhead_height = unit(12, "mm"), 
                    arrowhead_width = unit(6, "mm"), 
                    arrow_body_height = grid::unit(6, "mm")) +
    geom_text(aes(x = c(5246000,5290000),
                  y = "chr16",
                  label = c("AhDODAα1","AhCYP76AD2")),
              size = 11,
              nudge_y = 2.5) +
    coord_cartesian(ylim = c(0,4)) +
    scale_x_continuous(breaks = c(5250000, 5275000)) +
    theme(legend.position = "none",
          plot.margin = unit(c(0, 2, 0.5, 2), "cm"),
          axis.line = element_line(linewidth = 2),
          axis.ticks = element_line(linewidth = 1.5),
          axis.text.x = element_text(size=20),
          panel.grid.major.y = ggplot2::element_line(colour = "grey", 
                                                     linewidth = 1),
          #axis.title.x = element_text(size=40),
          axis.ticks.length = unit(.25, "cm"),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank()) +
    scale_fill_manual(values = c("chocolate2","cyan3",'red'))
  
  # combine plots:
  out <- plot_grid(p2,p1,
                   nrow = 2,
                   rel_heights = c(0.3,0.7))
  
  out
}

# plot all results
plot_AM_00332_leaf_green_red <- plotGresults1(AM_00332_leaf_green_red,
                                             betalain_genes = betalain_genes)
plot_AM_00331_flower_red_green <- plotGresults(AM_00331_flower_red_green,
                                               betalain_genes = betalain_genes)




plotleaf16 <- plotGqtl1(AM_00332_leaf_green_red,
                        genes = betalain_genes,
                        chr = 16)

plotflower16 <- plotGqtl(AM_00331_flower_red_green,
                         genes = betalain_genes,
                         chr = 16)
plotflower16


# use the cowplot package
cowplot_flower <- plot_grid(plot_AM_00331_flower_red_green, plotflower16,
                          labels = c("A", "B"),
                          nrow = 1,
                          align = "h",
                          rel_widths = c(0.7, 0.3),
                          label_size = 34)

ggsave(filename = "paper_grid_flower.png",
       plot = cowplot_flower,
       dpi = 400,
       width = 25,
       height = 5,
       bg = "white")

cowplot_leaf <- plot_grid(plot_AM_00332_leaf_green_red, plotleaf16,
                          labels = c("A", "B"),
                          nrow = 1,
                          align = "h",
                          rel_widths = c(0.7, 0.3),
                          label_size = 30)

ggsave(filename = "paper_grid_leaf.png",
       plot = cowplot_leaf,
       dpi = 400,
       width = 25,
       height = 5,
       bg = "white")


# combine with other plots:
pathway_plot <- ggdraw() + draw_image("plots/betalain_pathway_expression.png")

alignment_plot <- ggdraw() + draw_image("plots/AmMYB2_figure/S6_betalain_myb_R3.png")

MYB_plot <- plot_grid(cowplot_flower, pathway_plot, alignment_plot,
                       nrow = 3,
                       align = "v",
                      rel_heights = c(0.25, 0.75, 0.35),
                      labels = c("", "C", "D"),
                      label_size = 34)

ggsave(filename = "plots/paper_myb_combined_alignment.png",
       plot = MYB_plot,
       dpi = 400,
       width = 25,
       height = 25,
       bg = "white")







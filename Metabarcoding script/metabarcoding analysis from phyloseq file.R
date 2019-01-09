######Guillaume Nguyen Script phyloseq Visualisation version 1.0 (ready to use)
##### Filteau lab 3.12.18
#### guillaume.nguyen.quang@gmail.com for suggestion troubleshooting
## https://github.com/guillaumeQnguyen/Script for associated scripts (metagenomics...)
#

## libraries ####
library(reshape2)
library(rPackedBar)
library(phyloseq)
library(ggtree)
library(phangorn)
library(DECIPHER) 
library(Rmisc)
library(treemap)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)  
library(ggrepel)
library(colorspace)
library(scales)
library(PMCMR)
library (lmtest)
library(edgeR)
library(dendextend)
library(WGCNA)
# custom fonction ####
source("~/OneDrive - Universite Laval (1)/metagenomic/analyse Merge 16S ITS/Metagenomic function.R")






# Initial data ######
        setwd( "~/OneDrive - Universite Laval (1)/metagenomic/analyse Merge 16S ITS/")
        
        Phyloseq_file = readRDS( "2018-12-10filtrated 16S and ITS  merge.rds") #### read from phyloseq file
        subset =F ## for encoding test

## Extract data from phyloseq #####
        OTU.DF.main = data.frame(t(otu_table(Phyloseq_file))) ### extract OTU table
        table_sample = sample_data(Phyloseq_file) ## extract sample table
        table_sample[is.na(table_sample)] = ""
        taxa_table = data.frame(phyloseq::tax_table(Phyloseq_file),stringsAsFactors = F)
        taxa_table[is.na(taxa_table)] = ""
        taxa_table[is.na(taxa_table)] = ""
### Taxa formatage#####
                for (i in 1:nrow ( taxa_table)){
                if(taxa_table$Genus[i] == ""){
                  
                  taxa_table$Genus[i][taxa_table$Genus[i]==""] = taxa_table$Family[i]
                  taxa_table$Genus[i][taxa_table$Genus[i]==""] = taxa_table$Order[i]
                  taxa_table$Genus[i][taxa_table$Genus[i]==""] = taxa_table$Class[i]
                  taxa_table$Genus[i][taxa_table$Genus[i]==""] = taxa_table$Phylum[i]
                    taxa_table$Genus[i][taxa_table$Genus[i]==""] = taxa_table$Kingdom[i]
                } }
#### Subset application ######
       if (subset ==T ) {
  OTU.DF = OTU.DF.main[,1 :500 ]  ### make a small subset of data to test script
}else{OTU.DF = OTU.DF.main}

##### Melt for plot visualisation and merge with attribute ####
        OTU.DF$name = rownames( OTU.DF)
        M_OTU =  melt (OTU.DF, id.vars = "name" )
        table_sample_DF = data.frame(table_sample)
        M_OTU = merge (  M_OTU, table_sample_DF, by.x = "name", by.y = "Code.échantillon" )
        M_OTU = merge (  M_OTU, taxa_table, by.x = "variable", by.y = "OTU" )

###### Plot abundance fungi and Bacteria #####
      # Sum abundance Df for ordering #
          Sum_abundance =  aggregate(value ~name, M_OTU,FUN=sum)
          Sum_abundance = Sum_abundance[order( Sum_abundance$value,decreasing = F),]
          Sum_abundance$order = 1: nrow(Sum_abundance)
          Sum_abundance$name = factor ( Sum_abundance$name , Sum_abundance$name [Sum_abundance$order])
          M_OTU$name = factor (M_OTU$name)
          M_OTU$name2 =    factor(M_OTU$name, levels=levels(Sum_abundance$name)) ### ordering by max to min abundance sum
          
      ## Use custom function for ploting in log. Main Dataframe for next function. ##
          M_OTU_TP_b = proportion_DF(data_frame_OTU = M_OTU,Kingdom = "Bacteria") ## bacteria 
          M_OTU_TP_f  = proportion_DF(data_frame_OTU = M_OTU,Kingdom = "k__Fungi") ## fungi 
          
      ### Plot, style and order can be change ###
      
        g2 = ggplot(M_OTU_TP_b, aes (x = name, y=  proportion, fill = Genus)) + coord_flip() + geom_bar (stat = "identity")+ theme_bw() +
          theme(axis.title.x = element_blank(), 
                axis.title.y = element_blank(), 
                axis.text.y = element_blank(), 
                axis.ticks.y = element_blank(), 
                plot.margin = unit(c(1,-1,1,0), "mm") , legend.position = "right")  + ggtitle("Bacteria")
        
        g1 = ggplot(M_OTU_TP_f, aes (x = name, y=  proportion, fill = Genus)) + coord_flip() + geom_bar (stat = "identity")+ theme_bw() +
          theme(axis.title.x = element_blank(), 
                axis.title.y = element_blank(), 
                axis.text.y = element_blank(), 
                axis.ticks.y = element_blank(), 
                plot.margin = unit(c(1,-1,1,0), "mm") , legend.position = "left")  + ggtitle("Fungi") +  scale_y_reverse()

        M_OTU_T_plot = data.frame (name = unique (M_OTU[,2]))
        M_OTU_T_plot$name =  factor(M_OTU_T_plot$name , levels=levels(Sum_abundance$name))
        g.mid<-ggplot(M_OTU_T_plot,aes(x=1,y= name))+geom_text(aes(label=name))+
          geom_segment(aes(x=0.94,xend=0.96,yend=name))+
          geom_segment(aes(x=1.04,xend=1.065,yend=name))+
          ggtitle("")+
          ylab(NULL)+
          scale_x_continuous(expand=c(0,0),limits=c(0.94,1.065))+
          theme(axis.title=element_blank(),
                panel.grid=element_blank(),
                axis.text.y=element_blank(),
                axis.ticks.y=element_blank(),
                panel.background=element_blank(),
                axis.text.x=element_text(color=NA),
                axis.ticks.x=element_line(color=NA),
                plot.margin = unit(c(1,-1,1,-1), "mm"))
      
      #### Plot phyloseq basic merge ###
      gg1 <- ggplot_gtable(ggplot_build(g1))
      gg2 <- ggplot_gtable(ggplot_build(g2))
      gg.mid <- ggplot_gtable(ggplot_build(g.mid))
      grid.arrange(gg1,gg.mid,gg2,ncol=3,widths=c(4/9,1/9,4/9))
      
      # end plot abundance #

      
####### Plot prevalance #######
    prevalence_plot_phyloseq(taxa = "Family")
######## Species repartition ####
    # Plot abundance by species use custom function
      ## custom ccolor 
      my.cols <- brewer.pal(9, "Pastel1")
      my.cols[9] = "#9dd3a8"
      my.cols_Add = c ( "#c3a0c5", "#d39d9d",  "#9dd3c7", "#d0d39d","#a7d39d","#bd9dd3","#375587" , "#378738")
      my.cols = c ( my.cols, my.cols_Add)
      
      ### Bacteria pie
      
 
              pie_bacteria_sp_abundance =   plot_sp_abundance_pie(data_frame_OTU = M_OTU_TP_b, 
                                                        title = "Bacteria genus abundance",by_taxonomy = "Genus",
                                                         legend_position = "right" )
              pie_bacteria_sp_prevalence =   plot_sp_prevalence_pie(data_frame_OTU = M_OTU_TP_b, 
                                                                  title = "Bacteria genus prevalence",by_taxonomy = "Genus",
                                                                  legend_position = "left" )
              grid.arrange(pie_bacteria_sp_abundance,pie_bacteria_sp_prevalence,ncol=2,widths=c(1/2,1/2))
     
      ### Fungi pie        
              pie_fungi_sp_abundance =   plot_sp_abundance_pie(data_frame_OTU = M_OTU_TP_f, 
                                                        title = "Fungi genus abundance",by_taxonomy = "Genus",
                                                         legend_position = "right" )
              
              pie_fungi_sp_prevalence =   plot_sp_prevalence_pie(data_frame_OTU = M_OTU_TP_f, 
                                                                    title = "Fungi genus prevalence",by_taxonomy = "Genus",
                                                                    legend_position = "left" )
              grid.arrange(pie_fungi_sp_abundance,pie_fungi_sp_prevalence,ncol=2,widths=c(1/2,1/2))
              
######### Merge taxa by phylogeny########
          # OTU merge by similitude / identity (paramters can be changed)
          
              #bacteria (default parameter : similitude  = 4/465 , species with less than 10 OTU are not merge)
              group_consens_seq_b_c = NULL
              for ( b in unique (M_OTU_TP_b$Genus) ) {
                group_consens_seq =  plot_tree_sp(data.frame_OTU = M_OTU_TP_b,species = b,minimun_merge = 10)
                group_consens_seq_b_c = rbind (group_consens_seq_b_c ,group_consens_seq)
              }
              #fungi (default parameter : similitude  = 4/465 , species with less than 10 OTU are not merge)
              group_consens_seq_f_c = NULL
              for ( b in unique (M_OTU_TP_f$Genus) ) {
                group_consens_seq =  plot_tree_sp(data.frame_OTU = M_OTU_TP_f,species = b,minimun_merge =  10)
                group_consens_seq_f_c = rbind (group_consens_seq_f_c ,group_consens_seq)
              }
              Consensus_table_OTU =  rbind (group_consens_seq_f_c ,group_consens_seq_b_c)

              write.table (Consensus_table_OTU, "211218 table abundance metagenomic consensus.txt", quote =F, row.names=F,sep = "\t")
              
########## Phyloseq_remake_consensus ##### 
              cluster_OTU = dcast (Consensus  ~ name2  ,    data =   Consensus_table_OTU , value.var = "value"  ) #create matrix with mege OTU
              rownames(cluster_OTU) = cluster_OTU$Consensus
              cluster_OTU = cluster_OTU[,-1]
              cluster_OTU_value =  cluster_OTU
              cluster_OTU_phy = otu_table(data.frame(cluster_OTU),taxa_are_rows = T)
              taxa_consensus = unique (Consensus_table_OTU[,c (19:22 )] )
              rownames(taxa_consensus) = taxa_consensus$Consensus
              taxa_consensus_phy  = tax_table (as.matrix(taxa_consensus))
              sample_table_consenus  = unique (Consensus_table_OTU[,c (2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17 )] ) ## depend on the tablesample
              rownames(sample_table_consenus) = sample_table_consenus$name2
              sample_table_consenus_phy = sample_data (sample_table_consenus)
              consensus_phyloseq = merge_phyloseq(cluster_OTU_phy, taxa_consensus_phy,sample_table_consenus_phy )
              
           
########### Cluster repartition by WGCNA (see tutorial for more informations) #####
              binary=T  ### set for power determination if working with prevalence (T) or abundance(F) 
              if (binary == T){
                  Value_WGCNA = t(cluster_OTU_value)
                  Value_WGCNA[Value_WGCNA>1] = 1
              }else {  Value_WGCNA =  t(cluster_OTU_value) }
          #Check if sample are good to process 
            good= goodSamplesGenes(Value_WGCNA,verbose = 3) ### check is missing / error value for analysis
          ##dose the data is good to use
            good$allOK # if true you can continue
          ###distance tree  ### create dendro with attribution
            hclust(dist(Value_WGCNA)) #### make tree
          #### TOM script core
            powers = c ( c(1:10), seq(from =11, to = 30, by =2))
            sft = pickSoftThreshold(Value_WGCNA, powerVector = powers, networkType = "signed")
          ##### This is showing you the power (soft thresholding value), the r2 for the scale independence for 
          #each particular power (we shoot for an r2 higher than 0.8), the mean number of connections each node #
          #has at each power (mean.k), the median number of connections/node (median.k), and the maximum number 
          #of connections (max.k).
          ###### plot1
          cexl =0.9
          plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
               xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
               main = paste("Scale independence"));
          abline(h=0.8, col = "red")
          text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
               labels=powers,cex=cexl,col="red");
          
          ####### plot2
          plot(sft$fitIndices[,1], sft$fitIndices[,5],
               xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
               main = paste("Mean connectivity"))
          text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cexl,col="red")
          
    
############ Run wgcna on data ( abundance and prevalence, set your power threshold )#######
  # my data determinated with previous code  :  binary = 10, abundance = 19
          #abundance Group
            Value_WGCNA = t(cluster_OTU_value)
            DF_cluster_abun =  wgcna_home_function(phyloseq = consensus_phyloseq,OTU_table = Value_WGCNA,
                                       sample_attribution = sample_table_consenus,power_select =  10 ,set_minModuleSize = 3)
          ## prevalence Group
            Value_WGCNA[Value_WGCNA>1] = 1
            DF_cluster_preval =  wgcna_home_function(phyloseq = consensus_phyloseq,OTU_table = Value_WGCNA,
                                               sample_attribution = sample_table_consenus,power_select =  10 ,set_minModuleSize = 3)
          ### named col for merge
            colnames( DF_cluster_preval)[6] = "moduleColor_preval"
            colnames( DF_cluster_abun)[6] = "moduleColor_abun"
          #### Merging two group abundance and prevalence
            DF_cluster_merge = merge (DF_cluster_preval, DF_cluster_abun, by= "Taxon" )
          #####plot for difference OTU group between two method
             ggplot(DF_cluster_merge, aes( x = factor(moduleColor_preval),y = factor(moduleColor_abun)))+ geom_point(alpha = 0.2)+ 
             facet_wrap(~moduleColor_abun)
             
############# Plot by species #######
           # extract OTU by cluster
           DF_cluster_all = cluster_plot(data.frame = DF_cluster_preval, by.taxonomy = "Genus",
                                         DF_OTU_phyloseq = Value_WGCNA, sample_table= sample_table_consenus)
           
           
           ggplot (DF_cluster_all, aes (taxa, sub_group))+ geom_bar(stat= "identity")+ 
             facet_wrap(~ moduleColor_preval,scales = "free")
       
       
       
############## Plot diversity########

       plot_richness(Phyloseq_file, measures=c("Observed","Chao1"))
      
       plot_richness(Phyloseq_file, x = "Classement.2", color="Classement.2", measures=c("Chao1", "Shannon"))
      
###############  PLot ordination #####
      # define ordination
       out.pcoa.logt <- ordinate(Phyloseq_file, method = "PCoA", distance = "bray") ####PCoA or NMDS
       evals <- out.pcoa.logt$values$Eigenvalues
       
      ## plot sample
         plot_ordination(Phyloseq_file, out.pcoa.logt, 
                         color = "Producteur",label = "Code.échantillon") + labs(col = "Slash pile number") + 
         coord_fixed(sqrt(evals[2] / evals[1]))      
         plot_ordination(Phyloseq_file, out.pcoa.logt, type = "species", color = "Family") 
         p1 + facet_wrap(~Phylum, 5)
         
      ### plot species
        out.pcoa.logt <- ordinate(Phyloseq_file, "NMDS", "bray") ####PCoA or NMDS
        p1 = plot_ordination(Phyloseq_file, out.pcoa.logt, type = "species", color = "Family") 
        p1 + facet_wrap(~Phylum, 5)
     
      #### based on new method (plot)
      
   
       m = as(otu_table(Phyloseq_file), "matrix")
       # Add one to protect against overflow, log(0) issues.
       m = m + 1
       # Define gene annotations (`genes`) as tax_table
       taxonomy = tax_table(Phyloseq_file, errorIfNULL=FALSE)
       if( !is.null(taxonomy) ){
         taxonomy = data.frame(as(taxonomy, "matrix"))
       } 
       # Now turn into a DGEList
       d = DGEList(counts=m, genes=taxonomy, remove.zeros = TRUE)
       
       # Calculate the normalization factors
       z = calcNormFactors(d, method="RLE")
       # Check for division by zero inside `calcNormFactors`
       if( !all(is.finite(z$samples$norm.factors)) ){
         stop("Something wrong with edgeR::calcNormFactors on this data,
              non-finite $norm.factors, consider changing `method` argument")
       }
       
       plotMDS(z, col = as.numeric(factor(sample_data(Phyloseq_file)$Producteur)), labels = sample_names(Phyloseq_file))
  
################ End






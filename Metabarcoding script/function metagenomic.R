### function recall for metabarcoding data analysis


## function associated with data visualisation
proportion_DF  = function(data_frame_OTU =  M_OTU, Kingdom = "Bacteria"){
  
  M_OTU_T_b = data_frame_OTU[data_frame_OTU$Kingdom == Kingdom,]
  
  ###### bactria abbundance
  Sum_abundance_b =  aggregate(value ~name, M_OTU_T_b,FUN=sum)
  Sum_abundance_b = Sum_abundance_b[order( Sum_abundance_b$value,decreasing = F),]
  Sum_abundance_b$order = 1: nrow(Sum_abundance_b)
  Sum_abundance_b$name = factor ( Sum_abundance_b$name , Sum_abundance_b$name [Sum_abundance_b$order])
  
  #####axis order for plot
  M_OTU_T_b_s =  aggregate(value ~Genus, M_OTU_T_b,FUN=sum)
  M_OTU_T_b_s = M_OTU_T_b_s[order( M_OTU_T_b_s$value,decreasing = F),]
  M_OTU_T_b_s$order = 1: nrow(M_OTU_T_b_s)
  M_OTU_T_b_s$Genus = factor ( M_OTU_T_b_s$Genus , M_OTU_T_b_s$Genus [M_OTU_T_b_s$order])
  M_OTU_T_b$Genus =  factor(M_OTU_T_b$Genus, levels=levels(unique(M_OTU_T_b_s$Genus)))
  
  ####### log proportion for plot
  
  M_OTU_TP = M_OTU_T_b
  M_OTU_TP$proportion = 0
  
  for  (z in 1 : nrow (Sum_abundance_b )) {
    S1_sum = Sum_abundance_b[z,]
    S1_sum$log = log10 ( S1_sum$value)
    
    a = M_OTU_TP[ M_OTU_TP$name == as.character (S1_sum$name), ]
    
    ## logS = log t x abun / abun total
    
    M_OTU_TP$proportion[ M_OTU_TP$name == as.character (S1_sum$name)] =  S1_sum$log * a$value / S1_sum$value
  }
  
  M_OTU_TP$Genus =  factor(M_OTU_TP$Genus, levels=levels(M_OTU_T_b_s$Genus))
  M_OTU_TP$name = factor ( M_OTU_TP$name , Sum_abundance_b$name [Sum_abundance_b$order])
  
  return(M_OTU_TP)
}
plot_sp_abundance_pie =   function ( data_frame_OTU, title = "OTU sp abundance",legend_position = "left" ,by_taxonomy = "Genus"   ){
  taxonomy = which( colnames(data_frame_OTU) == by_taxonomy)
  Sum_sp_df = aggregate(value ~data_frame_OTU[,taxonomy], data_frame_OTU,FUN=sum)
  colnames(Sum_sp_df)[1]= "taxonomy"
  
  pie = ggplot(Sum_sp_df, aes ( x = "", y=  value, fill = taxonomy)) + 
    geom_bar (stat = "identity", )+ theme_bw() + coord_flip() +
    theme(axis.title.x = element_blank(), 
          axis.title.y = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank(), 
          plot.margin = unit(c(1,-1,1,0), "mm") , legend.position = legend_position)  + ggtitle(title) +labs(fill = by_taxonomy )
  pie = pie + coord_polar("y", start=0) + 
    geom_label_repel(size = 3,force = 3,aes(label = paste0(scientific(value,digits = 0))), position = position_stack(vjust = 0.5))
  +   scale_fill_manual(values = my.cols)
  #scale_fill_manual (values=colours)
  return ( pie)
}
plot_sp_prevalence_pie =   function ( data_frame_OTU, title = "OTU sp abundance",legend_position = "left" ,by_taxonomy = "Genus"   ){
  taxonomy = which( colnames(data_frame_OTU) == by_taxonomy)
  data_frame_OTU = data_frame_OTU[data_frame_OTU$name == data_frame_OTU$name[1], ]
  data_frame_OTU$prevalence = 1
  Sum_sp_df = aggregate(prevalence ~data_frame_OTU[,taxonomy], data_frame_OTU,FUN=sum)
  colnames(Sum_sp_df)[1]= "taxonomy"
  
  pie = ggplot(Sum_sp_df, aes ( x = "", y=  prevalence, fill = taxonomy)) + 
    geom_bar (stat = "identity", )+ theme_bw() + coord_flip() +
    theme(axis.title.x = element_blank(), 
          axis.title.y = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank(), 
          plot.margin = unit(c(1,-1,1,0), "mm") , legend.position = legend_position)  + ggtitle(title) +labs(fill = by_taxonomy )
  pie = pie + coord_polar("y", start=0) + 
    geom_label_repel(size = 3,force = 3,aes(label = prevalence), position = position_stack(vjust = 0.5))+
    scale_fill_manual(values = my.cols)
  #scale_fill_manual (values=colours)
  return ( pie)
}
plot_sp_abundance_pie_no_col =   function ( data_frame_OTU, title = "OTU sp abundance",legend_position = "left" ,by_taxonomy = "Genus"   ){
  taxonomy = which( colnames(data_frame_OTU) == by_taxonomy)
  Sum_sp_df = aggregate(value ~data_frame_OTU[,taxonomy], data_frame_OTU,FUN=sum)
  colnames(Sum_sp_df)[1]= "taxonomy"
  
  pie = ggplot(Sum_sp_df, aes ( x = "", y=  value, fill = taxonomy)) + 
    geom_bar (stat = "identity", )+ theme_bw() + coord_flip() +
    theme(axis.title.x = element_blank(), 
          axis.title.y = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank(), 
          plot.margin = unit(c(1,-1,1,0), "mm") , legend.position = legend_position)  + ggtitle(title) +labs(fill = by_taxonomy )
  pie = pie + coord_polar("y", start=0) 
  return ( pie)
}
plot_sp_prevalence_pie_no_col =   function ( data_frame_OTU, title = "OTU sp abundance",legend_position = "left" ,by_taxonomy = "Genus"   ){
  taxonomy = which( colnames(data_frame_OTU) == by_taxonomy)
  data_frame_OTU = data_frame_OTU[data_frame_OTU$name == data_frame_OTU$name[1], ]
  data_frame_OTU$prevalence = 1
  Sum_sp_df = aggregate(prevalence ~data_frame_OTU[,taxonomy], data_frame_OTU,FUN=sum)
  colnames(Sum_sp_df)[1]= "taxonomy"
  
  pie = ggplot(Sum_sp_df, aes ( x = "", y=  prevalence, fill = taxonomy)) + 
    geom_bar (stat = "identity", )+ theme_bw() + coord_flip() +
    theme(axis.title.x = element_blank(), 
          axis.title.y = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank(), 
          plot.margin = unit(c(1,-1,1,0), "mm") , legend.position = legend_position)  + ggtitle(title) +labs(fill = by_taxonomy )
  pie = pie + coord_polar("y", start=0)  + 
    geom_text(size = 1,aes(label = paste0(scientific(prevalence,digits = 0))))
  return (list ( Sum_sp_df, pie))
}


prevalence_plot_phyloseq  =    function ( Phyloseq = Phyloseq_file, taxa = "Phylum", preval.limit = 0.05 ){
  
  prevelancedf = apply(X = otu_table(Phyloseq),
                       MARGIN = 1,
                       FUN = function(x){sum(x > 0)})
  # Add taxonomy and total read counts to this data.frame
  prevelancedf = data.frame(Prevalence = prevelancedf,
                            TotalAbundance = taxa_sums(Phyloseq),
                            tax_table(Phyloseq))
  
  taxa_col = which ( colnames(prevelancedf) == taxa  )
  
  prevelancedf$select_taxa = prevelancedf[,taxa_col]
  
  prevelancedf1 = subset(prevelancedf, select_taxa %in% get_taxa_unique(Phyloseq, taxonomic.rank = taxa))
  
  
  plot =  ggplot(prevelancedf1, aes(TotalAbundance, Prevalence / nsamples(Phyloseq),color=select_taxa)) +
    # Include a guess for parameter
    geom_hline(yintercept = preval.limit, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
    scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
    facet_wrap(~select_taxa) + theme(legend.position="none")
  return(plot)
}

heatmap.plot.cluster_taxa = function (  correlation_matrix , adjusted.pvalue.matrix, 
                                        method.cluster = "ward.D2", k.cluster = 3 ,  name.write = species, title = "correlation",species_plot=species  ){
  
  select_name_write = name.write
  dend_r <- correlation_matrix %>% 
    dist %>% hclust(method =method.cluster ) %>% as.dendrogram  %>%
    color_branches(k=k.cluster)
  dend_c <- t(correlation_matrix) %>% 
    dist %>% hclust(method = method.cluster) %>% as.dendrogram  
  ##############cluster assignation
  # some_col_func <- function(n) rev(colorspace::heat_hcl(n, c = c(80, 30), l = c(30, 90), power = c(1/5, 1.5)))
  some_col_func <- colorspace::diverge_hcl
  # But due to the way heatmap.2 works - we need to fix it to be in the 
  # order of the data!   
  #rownames(d2_t)= 1 : nrow(d2_t)
  d2_t = d2
  col_labels <- col_labels[order(order.dendrogram(dend_r))]
  a =  rownames(matrix)[order.dendrogram(dend_r)]
  
  a = data.frame( a = a , order = 1: length(a))
  d2_tt =   merge(d2_t,a, by.x = "Taxon", by.y = "a")
  d2_tt = d2_tt[order(d2_tt$order),]
  
  col_labels2 = as.character(d2_tt$moduleColor)
  
  a$a != d2_tt$Taxon
  
  
  pdf(file = paste0( Sys.Date(), "/alignment/",by.taxonomy,species,"correlation.pdf"))
  gplots::heatmap.2(correlation_matrix, 
                    main = title,
                    srtCol = 35,
                    Rowv = dend_r,
                    Colv = dend_c,
                    trace="row", hline = NA,vline = NA , tracecol = NA,labRow = FALSE, labCol = colnames(correlation_matrix),        
                    margins =c(6,3),      
                    #key.xlab = "no / yes",
                    denscol = "grey",
                    density.info = "density",
                    col = some_col_func,
                    RowSideColors = col_labels2)
  
  dev.off()
}  

tip_glom_saveid <- function(physeq, h = 0.2, hcfun = agnes, ...){
  
  ## Necessary libraries
  require(phyloseq)
  require(ape)
  require(cluster)
  require(dplyr)
  require(tibble)
  
  ## Original phyloseq object
  physeq0 <- physeq
  
  ## From phyloseq::tip_glom
  dd = as.dist(cophenetic.phylo(phy_tree(physeq)))
  psclust = cutree(as.hclust(hcfun(dd)), h = h)
  #psclust = cutree(as.hclust(hcfun(dd, ...)), h = h)
  ## Key showing which OTUs are in which groups
  taxGlomKey <- data.frame(key = psclust) %>% rownames_to_column
  
  ## Key comparing new otus to old otus
  nGrp = max(taxGlomKey$key) 
  
  grps = vector("list", nGrp)
  for (loc_key in 1:nGrp){
    grps[[loc_key]] = taxGlomKey %>% filter(key == loc_key) %>% pull(rowname)
  }
  
  ## Turn list of lists into list of concatenated strings.
  oldGroups <- sapply(grps,
                      function(x){
                        unlist(
                          paste(x, collapse = ', ')
                        )
                      }
  )
  
  ## From phyloseq::tip_glom - do the actual merging
  cliques = levels(factor(psclust))
  consensus_seq = data.frame( cluster = 0, Consensus = "",stringsAsFactors = F)
  physeq2= physeq
  
  if(nGrp == 1) {
    seq_cluster = names(psclust)[psclust ==  1]
    consensus_seq =  as.character(ConsensusSequence(DNAStringSet(seq_cluster)))
    
    physeq2_seq_M = data.frame (cluster = 1 , Consensus = consensus_seq, rowname = seq_cluster  )}else{
      
      for (i in cliques) {
        
        seq_cluster = names(psclust)[psclust ==  i]
        
        physeq2 = merge_taxa(physeq2, eqtaxa = names(psclust)[psclust ==  i])
        
        consensus_seq[i,1] = i
        
        consensus_seq[i,2] =  as.character(ConsensusSequence(DNAStringSet(seq_cluster)))
        
      }
      physeq2_seq = data.frame (physeq2[3])
      physeq2_seq$cluster = rownames(physeq2_seq)
      physeq2_seq = merge ( physeq2_seq,consensus_seq, by.x = "cluster" )
      physeq2_seq= physeq2_seq[,-2]
      
      physeq2_seq_M = merge (physeq2_seq, by.x= "cluster",  taxGlomKey, by.y = "key" )
      
      
    }
  
  
  
  
  return (physeq2_seq_M)
}

plot_sp_preval_abun =   function( data.frame_OTU = M_OTU_TP_b, by.taxonomy = "Genus", species = "Pseudomonas", critere = "Producteur",
                                  cluster_annot = DF_cluster, cutoff_tree = 4, print_browse =F) { 
  species_plot1 = species
  taxonomy = which( colnames( data.frame_OTU) == by.taxonomy  )
  df_species =  data.frame_OTU[data.frame_OTU[,taxonomy] == species, ]
  ### prevalence in all sample
  df_species$prevalence_value = as.numeric(as.character(df_species$value))
  df_species$prevalence_value [df_species$prevalence_value>1]= 1
  df_species_preva_SUM = aggregate(prevalence_value ~variable ,df_species, FUN = sum  )
  colnames (df_species_preva_SUM)[2] = "prevalence"
  #### abundance in all sample
  Abun_species= aggregate(value ~variable, df_species,FUN=sum)
  colnames (Abun_species)[2] = "SPabundance"
  
  df_species_merge = merge (df_species , df_species_preva_SUM, by= "variable" )
  df_species_merge = merge (df_species_merge , Abun_species, by= "variable" )
  df_species_merge_S = unique(df_species_merge[, c (1,20:25,29,30)])
  
  library(RColorBrewer)
  
  myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
  
  p= ggplot ( df_species_merge_S, aes(x = SPabundance , y = prevalence ))+ theme_bw()  + stat_binhex()
  
  myColor_scale_fill_sqrt <- scale_fill_gradientn(colours = myPalette(4), trans = "sqrt")
  p1 <- p + myColor_scale_fill_sqrt + stat_binhex(bin = 5)
  P1 = p1 + ggtitle(species)
  P1
  
  
  ##### total count
  
  if (critere != "critere"){ 
    data.frame_SP = df_species
    col_select = which (colnames(data.frame_SP) == critere)
    
    df_species_preva_SUM_Critere = aggregate(prevalence_value ~variable + data.frame_SP[,col_select]  ,data.frame_SP, FUN = sum  )
    colnames (df_species_preva_SUM_Critere)[2:3] =  c("Critere","prevalence")
    
    
    Abun_species_critere= aggregate(value ~variable + data.frame_SP[,col_select], data.frame_SP,FUN=sum)
    colnames (Abun_species_critere)[2:3] = c("Critere","SPabundance")
    
    df_species_merge_critere = merge (data.frame_SP , df_species_preva_SUM_Critere,by.x =c ( "variable",critere)  , by.y= c ( "variable","Critere") )
    df_species_merge_critere = merge (df_species_merge_critere , Abun_species_critere, by.x =c ( "variable",critere)  , by.y= c ( "variable","Critere") )
    df_species_merge_S_critere = unique(df_species_merge_critere[, c (1,2,19,20:25,29,30)])
    
    library(RColorBrewer)
    
    myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
    
    col_select2 = which (colnames(df_species_merge_S_critere) == critere)
    
    df_species_merge_S_critere$log = log10 (df_species_merge_S_critere$SPabundance)
    df_species_merge_S_critere$log[df_species_merge_S_critere$log == "-Inf"] = 0
    
    
    p= ggplot ( df_species_merge_S_critere, aes(x = log , y = prevalence ))+ theme_bw()  +  facet_wrap(~df_species_merge_S_critere[,col_select2],scales = "free")
    
    myColor_scale_fill_sqrt <- scale_fill_gradientn(colours = myPalette(4), trans = "sqrt", breaks=c(0,5,20,50,100,200), limits=c(0,200))
    p1 <- p + myColor_scale_fill_sqrt + stat_binhex(bin = 5)
    P2 =  p1 + ggtitle(species)
  }else{P2 = ""}
  
  
  
  
  
  
  return ( list(P1,P2))
}      

wgcna_home_function = function ( phyloseq = Phyloseq_file, power_select = 10, OTU_table =Value_WGCNA, 
                                 sample_attribution = sample_table_consenus, merge_cluster = 0.30, set_minModuleSize= 5 )
{   
  #define your own powers limits according to your plot , r2 above 0.8 and the lowest power
            adj = adjacency( OTU_table, power = power_select, type = "signed")
            ### define topological overlap matrix
            Tom = TOMsimilarity(adj, TOMType = "signed")
            invTOM = 1-Tom
            TaxaTree = hclust(as.dist(invTOM), method = "average");
  ## save plot tree       
  pdf(file = paste0("distribution",".pdf"))
            plot(TaxaTree, xlab="", sub="", main = paste0("Taxa clustering on TOM-based dissimilarity power",power_select)   ,
                 labels = FALSE, hang = 0.04)
  dev.off()
  
  ####### test stats sur parametres#####
  ###### quantitatif
  # nPhenotype=38
  #network dendrogram with GS and module colors beneath
  # Calculate gene significance for Length
  # basePVal = 0.05;
  # z = qnorm(1-basePVal)/sqrt(nPhenotype-3);
  #baseCor = tanh(z);
  ########### Qualitative vs control.
  OTU_table.stats = data.frame((OTU_table))
  
  DF.pval = data.frame(OTU = colnames(OTU_table.stats),Sbb_pval =0 , Classement_pval = 0, Producteur =0,year = 0  )
  
  for ( x in 1: ncol(OTU_table.stats)){
    teh = data.frame(name = rownames(OTU_table.stats), OTU = OTU_table.stats[,x])
    select_OTU = base::merge (teh,sample_attribution, by.x = "name" , by.y = "name2" )
    select_OTU = data.frame (select_OTU)
    #fm1 =  lm(OTU ~ Sbb , select_OTU)
    #fg = lrtest(fm1)[5]
    #DF.pval[x,2] = fg[2,1]
    
    class_df = select_OTU[select_OTU$Classement.2 != "",]
    
    if( length(levels (factor(class_df$OTU))) == 1   ){
      DF.pval[x,3] = 1
    }else{
      
      fg1  = chisq.test(factor(class_df$OTU), factor(class_df$Classement.2))
      DF.pval[x,3] =  fg1[3][1]
    }
    fg2  = chisq.test(factor(select_OTU$OTU), factor(select_OTU$Producteur))
    DF.pval[x,4] =  fg2[3][1]
    
    fg3  = chisq.test(factor(select_OTU$OTU), factor(select_OTU$annee))
    DF.pval[x,5] =  fg3[3][1]
    
    
  }
  ####
  DF.pval$Sbbfdr = p.adjust(DF.pval$Sbb_pval,method = "fdr")
  DF.pval$Sbb_label = 0
  DF.pval$Sbb_label[ DF.pval$Sbbfdr< 0.05] = 1
  SbbpvalColors1 = numbers2colors( DF.pval$Sbb_label, signed = F);
  colnames(SbbpvalColors1) = "Sbb_pval"
  
  ######
  DF.pval$clafdr = p.adjust(DF.pval$Classement_pval,method = "fdr")
  DF.pval$cla_label = 0
  DF.pval$cla_label[ DF.pval$clafdr< 0.05] = 1
  clapvalColors1 = numbers2colors( DF.pval$cla_label, signed = F);
  colnames(clapvalColors1) = "classement_pval"
  
  
  ######
  DF.pval$profdr = p.adjust(DF.pval$Producteur,method = "fdr")
  DF.pval$pro_label = 0
  DF.pval$pro_label[ DF.pval$profdr< 0.05] = 1
  propvalColors1 = numbers2colors( DF.pval$pro_label, signed = F);
  colnames(propvalColors1) = "producteur_pval"
  
  ##########
  DF.pval$yearfdr = p.adjust(DF.pval$year,method = "fdr")
  DF.pval$year_label = 0
  DF.pval$year_label[ DF.pval$yearfdr< 0.05] = 1
  yearpvalColors1 = numbers2colors( DF.pval$year_label, signed = F);
  colnames(yearpvalColors1) = "year_pval"
  
  
  
  
  
  
  
  
  cor = cor(OTU_table.stats,as.numeric( as.character(sample_attribution$Sbb)), method = "spearman", use = "p"); ## change Sbb by quantitative 
  
  
  traitGeneColors1 = numbers2colors(cor, signed = TRUE);
  colnames(traitGeneColors1) = "Sbb";
  
  OTU_table.stats$name = rownames(OTU_table.stats)
  Stats_OTU_DF= merge (OTU_table.stats,by.x = "name", y = data.frame(sample_attribution),"name2"  )
  
  Stats_OTU_DF = Stats_OTU_DF[Stats_OTU_DF$Classement.2 != "",]
  
  Stats_OTU_DF$Classement.2 = factor ( Stats_OTU_DF$Classement.2)
  Stats_OTU_DF$Classement.2 =factor ( Stats_OTU_DF$Classement.2 , 
                                      levels( Stats_OTU_DF$Classement.2 )[c(5,6,1,2,3,4)])
  
  ###########
  
  minModuleSize = set_minModuleSize
  
  
  dynamicMods = cutreeDynamic(dendro = TaxaTree, distM = invTOM,
                              deepSplit = 2, pamRespectsDendro = FALSE,
                              minClusterSize = minModuleSize);
  table(dynamicMods)
  dynamicColors = labels2colors(dynamicMods)
  table(dynamicColors)
  
  sizeGrWindow(8,6)
  print(plotDendroAndColors(TaxaTree, dynamicColors, "Dynamic Tree Cut",
                            dendroLabels = FALSE, hang = 0.03,
                            addGuide = TRUE, guideHang = 0.05,
                            main = "Taxa dendrogram and module colors"))
  ###### similarity between cluster
  MEList = moduleEigengenes(OTU_table, colors = dynamicColors)
  MEs = MEList$eigengenes
  invMEs = 1-cor(MEs)
  MEDissThres = merge_cluster
  METree = hclust(as.dist(invMEs), method = "average")
  sizeGrWindow(7, 6)
  
  
  pdf(file = paste0("cluster tree",".pdf"))
  
  plot(METree, main = "Clustering of module eigengenes",
       xlab = "", sub = "")
  abline(h=MEDissThres, col = "red")
  
  
  dev.off()
  ##### merge cluster depending on distance
  merge = mergeCloseModules(OTU_table, dynamicColors, cutHeight = MEDissThres, verbose = 3)
  mergedColors = merge$colors;
  mergedMEs = merge$newMEs;
  sizeGrWindow(12, 9)
  
  pdf(file = paste0("cluster and tree",".pdf"))
  
  
  
  plotDendroAndColors(TaxaTree, cbind(dynamicColors, mergedColors,traitGeneColors1,SbbpvalColors1, clapvalColors1 , propvalColors1 , yearpvalColors1),
                      c("Dynamic Tree Cut", "Merged dynamic", "Sbb","Sbbpval", "Classement", "producteur","year"),
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  
  
  
  dev.off()
  
  
  ### retake nrege info for cluster 
  moduleColors = mergedColors
  colorOrder = c("grey", standardColors(50));
  moduleLabels = match(moduleColors, colorOrder)-1;
  MEs = mergedMEs;
  
  
  
  #################correlation with attribute
  nTaxa = ncol(OTU_table)
  nSamples = nrow(OTU_table)
  
  ###calculating
  
  MEs0 = moduleEigengenes(OTU_table, moduleColors)$eigengenes
  MEs = orderMEs(MEs0)
  moduleTraitCor = cor(MEs, sample_attribution[,c(2,1)], use = "p")
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
  
  sizeGrWindow(10,6)
  
  textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                     signif(moduleTraitPvalue, 1), ")", sep = "");
  dim(textMatrix) = dim(moduleTraitCor)
  par(mar = c(6, 8.5, 3, 3));
  
  pdf(file = paste0("attribute and cluster",".pdf"))
  
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = names(sample_attribution[,c(2,1)]),
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = FALSE,
                 colors = greenWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.5,
                 zlim = c(-1,1),
                 main = paste("Module-trait relationships"))
  dev.off()
  
  ####
  CR = as.data.frame(sample_attribution$Sbb);
  names(CR) = "CR"
  
  modNames = substring(names(MEs), 3)
  TaxaModuleMembership = as.data.frame(cor(OTU_table, MEs, use = "p"));
  MMPvalue = as.data.frame(corPvalueStudent(as.matrix(TaxaModuleMembership), nSamples));
  names(TaxaModuleMembership) = paste("MM", modNames, sep="");
  names(MMPvalue) = paste("p.MM", modNames, sep="");
  TaxaTraitSignificance = as.data.frame(cor(OTU_table, CR, use = "p"));
  GSPvalue = as.data.frame(corPvalueStudent(as.matrix(TaxaTraitSignificance), nSamples));
  names(TaxaTraitSignificance) = paste("GS.", names(CR), sep="");
  names(GSPvalue) = paste("p.GS.", names(CR), sep="");
  ####
  
  
  color_DF = data.frame ( moduleColors)
  probes = colnames(OTU_table)
  annot = data.frame(tax_table(phyloseq))
  probes2annot = match(probes, annot$Consensus)
  TaxaInfo0 = data.frame(Taxon = probes,
                         TaxaSymbol = annot$Consensus[probes2annot],
                         sub_group = annot$Sub.group[probes2annot],
                         taxa = annot$taxa[probes2annot],
                         cluster = annot$cluster[probes2annot],
                         moduleColor = moduleColors,
                         TaxaTraitSignificance,
                         GSPvalue)
  
  
  return(TaxaInfo0)
}

plot_tree_sp =   function( data.frame_OTU = M_OTU_TP_b,minimun_merge = 10, #have to be over 2
                           by.taxonomy = "Genus", species = "Rahnella"
                           , cutoff_tree = 4/465, print_browse =F, col.exp) { 
  title= paste ( "correlation",species)
  species_plot1 = species
  name.write = species
  taxonomy = which( colnames( data.frame_OTU) == by.taxonomy  )
  df_species =  data.frame_OTU[data.frame_OTU[,taxonomy] == species, ]
  
  ### prevalence in all sample
  df_species$prevalence_value = as.numeric(as.character(df_species$value))
  df_species$prevalence_value [df_species$prevalence_value>1]= 1
  
  
  seqs <- na.omit(unique(as.character(df_species$variable)))#getSequences(seqtabNoC)
  names(seqs) <- seqs # This propagates to the tip labels of the tree
  
  if( length( seqs) <= minimun_merge  ){
    
                      df_species$Consensus =  df_species$variable
                      df_species$Consensus= factor(df_species$Consensus)
                      df_species$cluster = "1"
                      df_species$Sub.group = as.numeric (df_species$Consensus)
                      df_species$taxa = df_species$Genus
                      
                      col_select =   which( colnames(df_species) %in% c(col.exp,"Consensus","prevalence_value","Sub.group", "taxa","cluster" ))  #c("value","cluster" , "Consensus"     ,"annee" ,"Sbb" ,    "Numéro.déchantillon..pairé.",
                                                                        #"Producteur","Classement.2","Station","Latitude","Longitude","Contenant","Date","Heure","Type","Brix..producteur.","Commentaires","Baril","name2","Sub.group","taxa", "prevalence_value"  ))
                      table_consensus = df_species[, col_select]
                      table_consensus$OTU_num = 1
                  }else{ 
                    alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA,verbose=FALSE)
                    
                    dir.create(path = paste0(  Sys.Date(), "/alignment" ))
                    write.table ( alignment, paste0( Sys.Date(), "/alignment/",by.taxonomy,species,".txt" ),sep = "\t",row.names = F,quote = F  )
                    writeXStringSet ( alignment, paste0( Sys.Date(), "/alignment/",by.taxonomy,species,"fname.fa"))

                    if (print_browse ==T){
                      BrowseSeqs(t, highlight=1)
                    }
                    d <- DistanceMatrix(alignment)
                    c <- IdClusters(d, method="complete", cutoff= cutoff_tree)
                    
                    names(alignment)<- seqs
                    phangAlign <- phyDat(as(alignment, "matrix"), type="DNA")
                    dm <- dist.ml(phangAlign)
                    treeNJ <- NJ(dm) # Note, tip order != sequence order
                    #fit = pml(treeNJ, data=phangAlign)
                    #fitGTR <- update(fit, k=4, inv=0.2)  ### default value
                    #### step take time
                    #fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE, 
                    #                  rearrangement = "stochastic", control = pml.control(trace = 0))
                    #plot(fitGTR)
                    #saveRDS ( fitGTR ,paste( "pseudomonas tree.pml",sep = ""))
                    #physeq = merge_phyloseq( fitGTR )
                    
                    physeq_plot = phyloseq(treeNJ)
                    
                    #physeq = phyloseq(GP.chl@otu_table,GP.chl@tax_table, GP.chl@sam_data, treeNJ)
                    
                    tree <- physeq_plot
                    p <- ggtree(tree, ladderize = F)
                    d1 <- data.frame(id=rownames(c), location= c$cluster ,stringsAsFactors = F) 
                    p1 <- p %<+% d1 + geom_tippoint(aes(color=factor(location)) , size = 2, position = "dodge",alpha = 0.7 )
                    P5 =p1 + ggtitle (paste (species, "by mismatch", cutoff_tree, "nucleotides"))
                    
                    
                    ############## plot taxo merge
                    
                    
                    test = tip_glom_saveid(physeq = physeq_plot,h = cutoff_tree)
                    
                    rownames(test)= test$rowname
                    
                    d2 <- data.frame(id=test$rowname, location= test$cluster )
                    p3 <- p %<+% d2 + geom_tippoint(aes(color=location) , size = 2, position = "dodge",alpha = 0.7 )
                    PP3 = p3 + ggtitle (paste (species, "OTU merge", cutoff_tree ))
                    
                    
                    correlation_Melt = test
                    
                    nam_file =  paste0( Sys.Date(), "/alignment/",by.taxonomy,species_plot1,"correlation.pdf")
                    pdf(file =nam_file)
                    print(PP3)
                    
                    f = na.omit ( unique (correlation_Melt$cluster))
                    for (h in f)
                                  {
                                    name_cluster = correlation_Melt$rowname[correlation_Melt$cluster ==h]
                                    df_species_cluster = df_species[df_species$variable %in%name_cluster,  ]
                                    
                                    matrix =na.omit( dcast(df_species_cluster ,   variable ~ name , value.var = "prevalence_value",sum))
                                    rownames(matrix)= matrix[,1]
                                    matrixabun =na.omit( dcast(df_species_cluster ,   variable ~ name , value.var = "value",sum))
                                    rownames(matrixabun)= matrixabun[,1]
                                    
                                    
                                    
                                    #######plot 
                                    if (nrow(matrix) <2){ }else{
                                      
                                      matrix = as.matrix(matrix[,-1])
                                      # heatmap.plot.cluster_taxa (matrix, species_plot = species_plot1)
                                      
                                      select_name_write = name.write
                                      dend_r <- matrix %>% 
                                        dist %>% hclust(method ="ward.D2" ) %>% as.dendrogram  %>%
                                        color_branches(k=1)
                                      dend_c <- t(matrix) %>% 
                                        dist %>% hclust(method = "ward.D2") %>% as.dendrogram  
                                      ##############cluster assignation
                                      # some_col_func <- function(n) rev(colorspace::heat_hcl(n, c = c(80, 30), l = c(30, 90), power = c(1/5, 1.5)))
                                      some_col_func <- colorspace::diverge_hcl
                                      
                                      
                                      gplots::heatmap.2(matrix, 
                                                        main = title,
                                                        srtCol = 35,
                                                        Rowv = dend_r,
                                                        Colv = dend_c,
                                                        trace="row", hline = NA,vline = NA , tracecol = NA,labRow = FALSE, labCol = colnames(matrix),        
                                                        margins =c(6,3),      
                                                        #key.xlab = "no / yes",
                                                        denscol = "grey",
                                                        density.info = "density",
                                                        col = some_col_func)
                                      
                                      matrixabun = as.matrix(matrixabun[,-1])
                                      # heatmap.plot.cluster_taxa (matrix, species_plot = species_plot1)
                                      
                                      select_name_write = name.write
                                      # dend_r <- matrixabun %>% 
                                      # dist %>% hclust(method ="ward.D2" ) %>% as.dendrogram  %>%
                                      #  color_branches(k=1)
                                      # dend_c <- t(matrixabun) %>% 
                                      #  dist %>% hclust(method = "ward.D2") %>% as.dendrogram 
                                      
                                      sum_taxa = apply(matrixabun,1, sum)
                                      sum_taxa_color =  numbers2colors(sum_taxa)
                                      
                                      ##############cluster assignation
                                      # some_col_func <- function(n) rev(colorspace::heat_hcl(n, c = c(80, 30), l = c(30, 90), power = c(1/5, 1.5)))
                                      some_col_func <- colorspace::sequential_hcl(length(matrix), c = c(100, 0), l = c(50, 90), power = 1)
                                      some_col_func = some_col_func[length(matrix):1]
                                      matrixabun_na = matrixabun
                                      matrixabun_na[matrixabun_na==0] = NA
                                      #### plot ####
                                      gplots::heatmap.2(matrixabun, 
                                                        main = title,
                                                        srtCol = 35,
                                                        Rowv = dend_r,
                                                        Colv = dend_c,
                                                        trace="row", hline = NA,vline = NA , tracecol = NA,labRow = FALSE, labCol = colnames(matrixabun),        
                                                        margins =c(6,3),      
                                                        #key.xlab = "no / yes",
                                                        denscol = "grey",
                                                        density.info = "density",
                                                        col = some_col_func,
                                                        RowSideColors = sum_taxa_color,  na.color = "white")
                                      
                                      
                        
                      }
    }
    dev.off()
    
    
    write.table ( unique(test$Consensus),paste0( Sys.Date(), "/alignment/",by.taxonomy,species,"consensus seq.txt"), quote = F, sep = "\t")
    
    
    
    ###### merge by consensus

    
    table_merge_consensus = merge ( test, by.x = "rowname" , df_species, by.y = "variable")
    table_merge_consensus$Sub.group = paste ( "Sub.group" , table_merge_consensus$cluster)
    table_merge_consensus[is.na(table_merge_consensus)] = ""
    
 
    table_consensus =  aggregate( value ~ cluster + Consensus + name2+Sub.group,table_merge_consensus ,  FUN = sum )
    
    table_merge_consensus$count = 1 
      
    table_merge_consensus_subg=  unique ( table_merge_consensus[,c (1,3,11,16)])
    table_merge_consensus_subg$OTU_num = 1 
    table_merge_consensus_subgS =  aggregate( OTU_num ~ Sub.group,table_merge_consensus_subg ,  FUN = sum )
    
    table_consensus = merge ( table_consensus, table_merge_consensus_subgS, by.x = c("Sub.group"), by.y = c("Sub.group"))
    
    # 
    # table_consensus = aggregate( value ~ 
    #                                cluster + Consensus     +annee +Sbb +    Numéro.déchantillon..pairé.+
    #                                Producteur+Classement.2+Station+Latitude+Longitude+Contenant+Date+Heure+Type+Brix..producteur.+Commentaires+Baril+name2+Sub.group                 
    #                              ,table_merge_consensus 
    #                              ,FUN = sum )
    # 
    
    
    table_consensus$prevalence_value = table_consensus$value
    table_consensus$prevalence_value[table_consensus$prevalence_value >=1] = 1
    
    nam_file2 =  paste0( Sys.Date(), "/alignment/",by.taxonomy,species_plot1,"merge OTU.pdf")
    pdf(file =nam_file2)
    
    
    matrix =na.omit( dcast(table_consensus ,   Sub.group ~ name2 , value.var = "prevalence_value",sum))
    rownames(matrix)= matrix[,1]
    matrixabun =na.omit( dcast(table_consensus ,   Sub.group ~ name2 , value.var = "value",sum))
    rownames(matrixabun)= matrixabun[,1]
    
    matrix = as.matrix(matrix[,-1])
    # heatmap.plot.cluster_taxa (matrix, species_plot = species_plot1)
    
    if (nrow(matrix == 1)){
      
      matrix_m =data.frame( t (matrix))
      matrix_m$name =  rownames(matrix_m)
      ggplot(matrix_m, aes( x = name, y = "seq" )) + geom_tile ( aes( fill =matrix_m[,1]))
      
    }else{
      select_name_write = name.write
      dend_r <- matrix %>% 
        dist %>% hclust(method ="ward.D2" ) %>% as.dendrogram  %>%
        color_branches(k=1)
      dend_c <- t(matrix) %>% 
        dist %>% hclust(method = "ward.D2") %>% as.dendrogram  
      ##############cluster assignation
      # some_col_func <- function(n) rev(colorspace::heat_hcl(n, c = c(80, 30), l = c(30, 90), power = c(1/5, 1.5)))
      some_col_func <- colorspace::diverge_hcl
      
      
      gplots::heatmap.2(matrix, 
                        main = title,
                        srtCol = 35,
                        Rowv = dend_r,
                        Colv = dend_c,
                        trace="row", hline = NA,vline = NA , tracecol = NA,labRow = FALSE, labCol = colnames(matrix),        
                        margins =c(6,3),      
                        #key.xlab = "no / yes",
                        denscol = "grey",
                        density.info = "density",
                        col = some_col_func)
      
      matrixabun = as.matrix(matrixabun[,-1])
      # heatmap.plot.cluster_taxa (matrix, species_plot = species_plot1)
      
      select_name_write = name.write
      # dend_r <- matrixabun %>% 
      # dist %>% hclust(method ="ward.D2" ) %>% as.dendrogram  %>%
      #  color_branches(k=1)
      # dend_c <- t(matrixabun) %>% 
      #  dist %>% hclust(method = "ward.D2") %>% as.dendrogram 
      
      sum_taxa = apply(matrixabun,1, sum)
      sum_taxa_color =  numbers2colors(sum_taxa)
      
      ##############cluster assignation
      # some_col_func <- function(n) rev(colorspace::heat_hcl(n, c = c(80, 30), l = c(30, 90), power = c(1/5, 1.5)))
      some_col_func <- colorspace::sequential_hcl(length(matrixabun), c = c(100, 0), l = c(50, 90), power =1)
      some_col_func = some_col_func[length(matrix):1]
      #some_col_func[1] = "#FFFFFF"
      #some_col_func[2] = "#000000"
      matrixabun_na = matrixabun
      matrixabun_na[matrixabun_na==0] = NA
      gplots::heatmap.2(matrixabun_na, 
                        main = title,
                        srtCol = 35,
                        Rowv = dend_r,
                        Colv = dend_c,
                        trace="row", hline = NA,vline = NA , tracecol = NA,labRow = FALSE, labCol = colnames(matrixabun),        
                        margins =c(6,3),      
                        #key.xlab = "no / yes",
                        denscol = "grey",
                        density.info = "density",
                        col = some_col_func,
                        RowSideColors = sum_taxa_color,
                        na.color = "white")
      
      
      
    }
    
    dev.off()
    
    table_consensus$taxa = species_plot1
   # table_consensus$OTU_num = length( seqs)
    
    
    
  } 
  return (table_consensus)
  }

heatmap.plot.cluster = function (  correlation_matrix , adjusted.pvalue.matrix, 
                                   method.cluster = "ward.D2", k.cluster = 3 ,  name.write = "sample", title = "correlation"  ){
  
  select_name_write = name.write
  dend_r <- correlation_matrix %>% 
    dist %>% hclust(method =method.cluster ) %>% as.dendrogram  %>%
    color_branches(k=k.cluster)
  dend_c <- t(correlation_matrix) %>% 
    dist %>% hclust(method = method.cluster) %>% as.dendrogram  %>%
    color_branches(k=k.cluster)
  ##############cluster assignation
  # some_col_func <- function(n) rev(colorspace::heat_hcl(n, c = c(80, 30), l = c(30, 90), power = c(1/5, 1.5)))
  some_col_func <- colorspace::diverge_hcl
  gplots::heatmap.2(correlation_matrix, 
                    main = title,
                    srtCol = 35,
                    Rowv = dend_r,
                    Colv = dend_c,
                    trace="row", hline = NA,vline = NA , tracecol = NA,labRow = FALSE, labCol = colnames(correlation_matrix),        
                    margins =c(6,3),      
                    #key.xlab = "no / yes",
                    denscol = "grey",
                    density.info = "density",
                    col = some_col_func )
  
}


cluster_plot = function (data.frame = DF_cluster_preval, DF_OTU_phyloseq = Value_WGCNA,
                         sample_table=sample_table_consenus, by.taxonomy = "Genus" ){
  
  
  
  
  #### format OTU df count
  OTU.DF_select_cluster = DF_OTU_phyloseq#[,name_cluster ]
  OTU.DF_select_cluster = data.frame (OTU.DF_select_cluster,stringsAsFactors = F)
  OTU.DF_select_cluster$name = rownames(OTU.DF_select_cluster)
  OTU.DF_select_cluster = melt(OTU.DF_select_cluster, id.var = 'name')
  OTU.DF_select_cluster$name = factor(OTU.DF_select_cluster$name)
  
  ####
  
  #### merge df
  OTU.DF_select_cluster =base::merge (x =OTU.DF_select_cluster, y = sample_table, by.x = "name" , by.y = "name2",all = T )
  OTU.DF_select_cluster_Msam_clust = merge (x=OTU.DF_select_cluster, y= data.frame , by.x = "variable", by.y = "Taxon")
  
  
  OTU.DF_select_cluster_Msam_clust =na.omit(OTU.DF_select_cluster_Msam_clust)
  OTU.DF_select_cluster_Msam_clust$Sbb =as.numeric(as.character (OTU.DF_select_cluster_Msam_clust$Sbb))
  OTU.DF_select_cluster_Msam_clust$Sbb2 = cut(OTU.DF_select_cluster_Msam_clust$Sbb, 
                                              breaks = c(0,25,50,100 ,300),
                                              labels = c("[0,25]","(25,50]","(50,100]","(100,300]" ),
                                              right = T)
  OTU.DF_select_cluster_Msam_clust$log = log10(OTU.DF_select_cluster_Msam_clust$value)
  OTU.DF_select_cluster_Msam_clust$log[OTU.DF_select_cluster_Msam_clust$log==-Inf] = 0
  
  print (ggplot (OTU.DF_select_cluster_Msam_clust, aes(x = Sbb , y = value, group = variable))+ 
           geom_line(alpha = 0.2,  stat = "smooth", method = "loess") +geom_point()+
           theme(legend.position="") + facet_wrap(~moduleColor_preval))
  
  
  print(ggplot(OTU.DF_select_cluster_Msam_clust, aes(x = name , fill = value, y = variable))+ geom_tile() + theme( axis.text.y=element_blank()) )
  nam_file =  paste0( Sys.Date(), "/alignment/",by.taxonomy,"Cluster correlation cluster.pdf")
  
  pdf(file =nam_file)
  for ( x in unique  (  data.frame$moduleColor)){ 
    select_cluster = data.frame[data.frame$moduleColor == x,]
    select_cluster.DF =  Value_WGCNA[,select_cluster$Taxon]
    # matrix = dcast(select_cluster.DF ,   variable ~ name , value.var = "value")
    # rownames(matrix)= matrix[,1]
    #matrix = as.matrix(matrix[,-1])
    heatmap.plot.cluster(correlation_matrix = select_cluster.DF,k.cluster = 1,name.write = "Test")
  }
  dev.off()
  
  return (OTU.DF_select_cluster_Msam_clust)
}


circleplot = function (data = Consensus_table_OTU, 
                       by.col.taxa = "taxa",by.col.value= "OTU_num"  ,by.col.group = "Sub.group",
                       type = "OTU", palette = "RdPu", police.range = c(5,10)  ){
  
  col_taxa = which ( colnames ( data) == by.col.taxa)
  col_value = which ( colnames ( data) == by.col.value)
  col_group = which ( colnames ( data) == by.col.group)
  
  data_u = unique ( data [, c (col_taxa, col_group, col_value  )])
  
  data_u$target = paste( data_u$taxa, data_u$Sub.group, sep = "_")
  
  colnames(data_u) = c ( "source", "0target","N","target")
  
  
  #data_u = data_u[data_u$source == "Pseudomonas",]
  
  # Generate the layout. This function return a dataframe with one line per bubble. 
  # It gives its center (x and y) and its radius, proportional of the value
  # We need a data frame giving a hierarchical structure. Let's consider the flare dataset:
  #edges=flare$edges
  edges = data.frame(from = as.character(data_u$source), to = as.character(data_u$target))
  edges = rbind( edges, data.frame(from = "taxa" ,  to =  unique(as.character(data_u$source)) ,stringsAsFactors = T))
  rownames(edges) = as.character(edges$to)
  # Usually we associate another dataset that give information about each node of the dataset:
  #vertices = flare$vertices
  vertices = data.frame(name =as.character(data_u$target), size=data_u$N, shortName =data_u$N)
  vertices = rbind( vertices, data.frame(name = unique(as.character(data_u$source)) ,  size = 1, shortName= unique(as.character(data_u$source))  ))
  vertices = rbind( vertices, data.frame(name = "taxa" ,  size = 1, shortName= "taxa" ))
  
  rownames(vertices) = as.character(vertices$name)
  
  # Then we have to make a 'graph' object using the igraph library:
  tree <- FromDataFrameNetwork(edges)
  mylevels=data.frame( name=tree$Get('name'), level=tree$Get("level") )
  vertices = vertices %>% left_join(., mylevels, by=c("name"="name"))
  
  vertices = vertices %>% mutate(new_label=ifelse(level==2, shortName, NA))
  
  mygraph <- graph_from_data_frame( edges, vertices=vertices )
  
  # Make the plot
  ggraph(mygraph, layout = 'circlepack', weight="size" ) + 
    geom_node_circle(aes(fill = depth)) +
    geom_node_text( aes(label=shortName, filter=leaf, fill=depth, size=size)) +
    geom_node_label( aes(label=new_label), size=3,alpha=0.7) +
    theme_void() + 
    theme(legend.position="FALSE") + 
    scale_fill_viridis() +
    scale_size_continuous(range =  police.range)  + scale_fill_distiller(palette = palette) 
  
}

bubbleplot_prevalence =  function (Df = Sum_sp_df,by.col.taxa = "taxonomy",by.col.value= "prevalence",  limit = 1, radius_size = 0.99,size = 1, police_size_min = 2, police_size_max= 5  ) {
  
  col_taxa = which ( colnames ( Df) == by.col.taxa)
  col_value = which ( colnames ( Df) == by.col.value)
  Df = Df[Df[,col_value] >limit,]
  Df[,col_taxa] = as.character( Df[,col_taxa])
  packing <- circleProgressiveLayout(log2 (Df[,col_value])+size, sizetype='area')
  packing$radius=0.99*packing$radius
  
  # We can add these packing information to the initial data frame
  data = cbind(Df, packing)
  
  # Check that radius is proportional to value. We don't want a linear relationship, since it is the AREA that must be proportionnal to the value
  plot(data$radius, data$value)
  
  # The next step is to go from one center + a radius to the coordinates of a circle that
  # is drawn by a multitude of straight lines.
  dat.gg <- circleLayoutVertices(packing, npoints=50)
  
  data$taxonomy= sub ("_", "\n",data$taxonomy )
  data$taxonomy= sub ("_", "\n",data$taxonomy )
  data$taxonomy= sub ("-", "\n",data$taxonomy )
  data$taxonomy= sub ("-", "\n",data$taxonomy )
  data$taxonomy= sub ("_", "\n",data$taxonomy )
  data$taxonomy= sub ("_", "\n",data$taxonomy )
  data$taxonomy= sub ("-", "\n",data$taxonomy )
  data$taxonomy= sub ("-", "\n",data$taxonomy )
  # Make the plot
  P = ggplot() + 
    
    # Make the bubbles
    geom_polygon(data = dat.gg, aes(x, y,  fill=as.factor(id)), colour = "black", alpha = 0.6) +
    
    # Add text in the center of each bubble + control its size
    geom_text(data = data, aes(x, y, label=paste (Df[,col_taxa],"\n", Df[,col_value]), size = log2(Df[,col_value]))) +
    scale_size_continuous(range = c(police_size_min,police_size_max)) +
    
    # General theme:
    theme_void() + 
    theme(legend.position="none") +
    coord_equal()
  
  return (P)
}

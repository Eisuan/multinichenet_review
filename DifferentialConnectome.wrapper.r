# Defensive checks for ensuring correct data processing 
# in function DifferentialConnectome.wrapper

test_equal_labels<-function(label1,label2){
  
  # in the wrapper, each label will have added _condition
  # remove _ and test if you are comparing the same cell labels
  # across different conditions
  label1_sub<-gsub("_(.*)", "",label1)
  label2_sub<-gsub("_(.*)", "",label2)
  
  if(label1_sub!=label2_sub){
    stop("Same cell types are not being correctly compared for DEA. Please check
         the code of DifferentialConnectome.wrapper function")
  }
  
}





# Define the Differential connectome wrapperbased on vignette

DifferentialConnectome.wrapper<-function(seurat_obj, metadata_cond_name, 
                                         test_cond, cntrl_cond,
                                         subset_obj=TRUE, db=NA, lognorm=TRUE){
  
  subset_ids=c(test_cond, cntrl_cond)
  
  # input validation
  if(subset_obj==TRUE){
    subset_oi_logical<-seurat_obj@meta.data[,metadata_cond_name]%in%subset_ids 
    subset_oi_cellid<-rownames(seurat_obj@meta.data[subset_oi_logical,])
    
    seurat_obj<-subset(seurat_obj, cells =subset_oi_cellid)
    }
  
  
  
  
  # Retriving genes in the default DB to perform filtering
  
  if(!is.data.frame(db)){                                                           # DEBUG!
     connectome.genes <- union(Connectome::ncomms8866_human$Ligand.ApprovedSymbol,
                            Connectome::ncomms8866_human$Receptor.ApprovedSymbol)
     genes <- connectome.genes[connectome.genes %in% rownames(seurat_obj)]
  }else{
     db.genes <- union(db[,1],db[,2])
     genes <- db.genes[db.genes %in% rownames(seurat_obj)] 
    
  }
  
  # Split the object by condition:
  seurat_subset_sce.list <- SplitObject(seurat_obj ,split.by = metadata_cond_name)
  
  
  # Normalize, Scale, and create Connectome:
  seurat_subset_sce.con.list <- list()
  for (i in 1:length(seurat_subset_sce.list)){
    
    
    if(lognorm==TRUE){
      seurat_subset_sce.list[[i]] <- NormalizeData(seurat_subset_sce.list[[i]])
    }
    seurat_subset_sce.list[[i]] <- ScaleData(seurat_subset_sce.list[[i]],
                                             features = rownames(seurat_subset_sce.list[[i]]))
    
    
    ## Create connectomes
    
    if(!is.data.frame(db)){
      seurat_subset_sce.con.list[[i]] <- CreateConnectome(seurat_subset_sce.list[[i]],
                                                        species = 'human',p.values = F)
    }else{
      seurat_subset_sce.con.list[[i]] <- CreateConnectome(seurat_subset_sce.list[[i]],
                                                          species = 'human',p.values = F,
                                                          LR.database="custom",
                                                          custom.list=db)
    }
  }
  names(seurat_subset_sce.con.list) <- names(seurat_subset_sce.list)
  
  # free up RAM space

  
  
  diff_E<-DifferentialConnectome(seurat_subset_sce.con.list[[cntrl_cond]], 
                         seurat_subset_sce.con.list[[test_cond]] )
  
  
  
  
  
  # Stash idents and make new identities which identify each as stimulated vs. control
  
  # subset seurat objs

  
  celltypes_e <- as.character(unique(Idents(seurat_obj)))
  celltypes_e.stim <- paste(celltypes_e, test_cond , sep = '_')
  celltypes_e.ctrl <- paste(celltypes_e, cntrl_cond, sep = '_')
  seurat_obj$celltype.condition <- paste(Idents(seurat_obj), 
                                            seurat_obj@meta.data[,metadata_cond_name], sep = "_")
  seurat_obj$celltype <- Idents(seurat_obj)
  Idents(seurat_obj) <- "celltype.condition"
  
  
  
  
  
  
  # Identify which ligands and receptors, for which cell populations, have an adjusted p-value < 0.05 based on a Wilcoxon rank test
  diff.p_e <- data.frame()
  for (i in 1:length(celltypes_e)){
    
    # test correct DEA between cell types of the same types among conditions
    test_equal_labels(celltypes_e.stim[i], celltypes_e.ctrl[i])
    
    temp <- FindMarkers(seurat_obj, 
                        ident.1 = celltypes_e.stim[i], 
                        ident.2 = celltypes_e.ctrl[i],
                        verbose = FALSE,
                        features = genes,
                        min.pct = 0,
                        logfc.threshold = 0)
    temp2 <- subset(temp, p_val_adj < 0.05)
    if (nrow(temp2)>0){
      temp3 <- data.frame(genes = rownames(temp2),cells = celltypes_e[i])
      diff.p_e <- rbind(diff.p_e, temp3)
    }
  }
  diff.p_e$cell.gene <- paste(diff.p_e$cells,diff.p_e$genes,sep = '.')
  
  # Filter differential connectome to only include significantly perturbed edges
  diff_E$source.ligand <- paste(diff_E$source,diff_E$ligand,sep = '.')
  diff_E$target.receptor <- paste(diff_E$target,diff_E$receptor,sep = '.')
  
  diff_E.subset <- subset(diff_E,source.ligand %in% diff.p_e$cell.gene & target.receptor %in% diff.p_e$cell.gene)
  
  diff_E.subset
  
}






build.connectome<-function(seurat_obj, db=NA, p.values = TRUE ){
   seurat_obj <- NormalizeData(seurat_obj)
   
   if(!is.data.frame(db)){                                                           # DEBUG!
     connectome.genes <- union(Connectome::ncomms8866_human$Ligand.ApprovedSymbol,
                               Connectome::ncomms8866_human$Receptor.ApprovedSymbol)
     genes <- connectome.genes[connectome.genes %in% rownames(seurat_obj)]
   }else{
     db.genes <- union(db[,1], db[,2])
     genes <- db.genes[db.genes %in% rownames(seurat_obj)] 
     
   }   
   seurat_obj <- ScaleData(seurat_obj,features = genes)


   
   
   if(!is.data.frame(db)){
     seurat_obj.con <- CreateConnectome(seurat_obj, species = 'human',p.values = p.values)
   }else{
     seurat_obj.con <- CreateConnectome(seurat_obj,species = 'human',p.values =p.values ,
                                                         LR.database="custom",
                                                         custom.list=db)
   }
   
   return(seurat_obj.con)
}





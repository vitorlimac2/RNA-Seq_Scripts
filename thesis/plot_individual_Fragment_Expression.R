  
### Medir a distancia da expressão de genes isolados contra um background aleatorio
  
  ## Le o gene alvo com sua distribuicao de fragmentos; 
  ## os genes de background com sua distribuicao de fragmentos;
  ## TPM normalizado e não normalizado
  
  ## Expressão do backgroun da replica 1 x replica 2
  ## Expressão do gene alvo da replica 1 x replica 2
  
  
  # 1 - Plot dos fragmentos do gene alvo contra background
  
  expression.comparison.plot <- function(target_genes, background_genes, raw_expression_table, norm_expression_table, replicate_names, target_replicate){
    
    if(is.null(target_replicates))
      target_replicates <- replicate_names
    
    # Plot expressão não normalizada, expressão normalizada e distribuição de fragmentos
    ncols <- 3
    nrows <- length(target_replicates)
    
      ## pra cada gene alvo, plota do fragmentos
    
    for (row in 1:nrow(target_genes)) {
      
      target_gene_id <- rownames(target_genes)[row]
      background_genes_id <- rownames(background_genes)
      
      ## plotar todas as expressões juntas
      par(mfrow=c(nrows,ncols))
      
      bg_raw <- log10(raw_expression_table[background_genes_id,]+1)
      bg_norm <- log10(norm_expression_table[background_genes_id,]+1)
      
      gene_raw <- log10(raw_expression_table[target_gene_id,]+1)
      gene_norm <- log10(norm_expression_table[target_gene_id,]+1)
     
      colnames(bg_raw) <- replicate_names
      colnames(bg_norm) <- replicate_names
      colnames(gene_raw) <- replicate_names
      colnames(gene_norm) <- replicate_names
      
      n <- length(target_replicates)
      
      for(i in 1:n){
        for(j in 1:length(replicate_names)){
          if(target_replicate[i] == replicate_names[j])
            next;
  
          plot(main = "Raw counts", 
               bg_raw[,target_replicates[i]],
               bg_raw[,replicate_names[j]], 
               col="black", 
               type = "p",
               xlab = paste(target_replicates[i],"(ln+1)"),
               ylab = paste(replicate_names[j],"(ln+1)"),
               pch=21, bg="lightgray")
          
          # adicionar expressao normalizado
          points(gene_raw[,target_replicates[i]],
                 gene_raw[,replicate_names[j]], col = "red", pch=16)
          
          points(gene_norm[,target_replicates[i]],
                 gene_norm[,replicate_names[j]], col = "blue", pch=16)
          
          
          legend("bottomright", 
                 legend=c(target_gene_id,"Background genes"), 
                 col=c("red", "black"), inset = 0.02,box.lty=0,
                 cex = 0.6)
          
          
          # Plot da expressao dos genes de background
          rep_name <- replicate_names
          plot(main = "Normalized counts", 
               bg_norm[,target_replicates[i]],
               bg_norm[,replicate_names[j]], 
               col="black", 
               type = "p",
               xlab = paste(target_replicates[i],"(ln+1)"),
               ylab = paste(replicate_names[j],"(ln+1)"),
               pch=21, bg="lightgray")
          
          points(gene_raw[,target_replicates[i]],
                 gene_raw[,replicate_names[j]], col = "red", pch=16)
          
          
          # adicionar o valor de expressao nao normalizado
          points(gene_norm[,target_replicates[i]],
                 gene_norm[,replicate_names[j]], col = "blue", pch=16)
          
  
          
          ## Do lado, colocar o plot da fragmentacao
          
          fragment.plot(target_gene_id,target_genes[target_gene_id,], background_genes)
        }
      }
    }
  }
  
  fragment.plot <- function(target_gene,target_dist, background_genes){
    
    ## inicializa vetor total de fragmentos
    total_dist <- strsplit(as.character(background_genes[1,1]),",")
    total_dist <- as.numeric(unlist(total_dist));
    
    ## armazena todos os valores de fragmentos
    for(row in 2:nrow(background_genes)){
      temp <- strsplit(as.character(background_genes[row,1]),",")
      total_dist <- c(total_dist,as.numeric(unlist(temp)));
    }
    
      geneid <- target_gene
      target_dist <- strsplit(as.character(target_dist),",")
      target_dist <- as.numeric(unlist(target_dist));
      
      my_breaks <- pretty(range(c(target_dist, total_dist)), n=10)
      
      h_target <- hist(target_dist, breaks=my_breaks, plot=FALSE)
      h_bg <- hist(total_dist, breaks=my_breaks, plot=FALSE)
      
      max_y <- max(h_target$counts*100/length(target_dist), h_bg$counts*100/length(total_dist))
      
      plot(main = "Fragment plot", h_target$mids, h_target$counts*100/length(target_dist), type="l", ylab = "%", xlab = "Read length (bp)", col="red", ylim=c(0,max_y)) 
      
      lines(h_bg$mids, h_bg$counts*100/length(total_dist), col="black")
      
      legend("topright", 
             legend=c(geneid,"Background genes"), 
             col=c("red", "black"), inset = 0.02,box.lty=0,
             cex = 0.5, lty = 1)
  }
  
  workdir <- "/home/vitor/Downloads"
  setwd(workdir)
  path_to_target_genes <- "ProR1.Target"
  #path_to_target_genes <- "ProR1.HighFragmentation.HighExpressed.100"
  #path_to_background_genes <- "ProR1.NormalFragmentationGenes"
  path_to_background_genes <- "ProR1.NormalFragmentation.Background.200"
  path_to_TPM_raw <- "TPM_raw_frac.csv"
  path_to_TPM_norm <- "TPM_normalized_frac.csv" 
  
  
  
  firstRepCol <- 3
  lastRepCol <- 15
  
  target_replicate <- "ProR1"
  
  replicate_names <- c("ProC1","ProC2","ProC3","ProC4","ProC5","ProC6","ProC7","ProC8","ProC9","ProR1","ProR2")
  
  target_replicates <- c("ProR1")
  
  target_genes <- read.table(path_to_target_genes, header = F, row.names = 1)
  background_genes <- read.table(path_to_background_genes, header = F, row.names = 1)
  raw_expression_table <- read.table(path_to_TPM_raw, header=FALSE, skip = 1, row.names = 1)[,firstRepCol:lastRepCol]
  norm_expression_table <- read.table(path_to_TPM_norm, header=FALSE, skip = 1, row.names = 1)[,firstRepCol:lastRepCol]
  
  expression.comparison.plot(target_genes,background_genes, 
                             raw_expression_table,
                             norm_expression_table,
                             replicate_names,
                             target_replicates)
  
  #################################################################################
  
  

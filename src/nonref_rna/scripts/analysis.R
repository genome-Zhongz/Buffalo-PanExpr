#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==3)
work <- normalizePath(args[1],mustWork=TRUE); out <- normalizePath(args[2],mustWork=TRUE)
if(!requireNamespace('jsonlite',quietly=TRUE))stop('Install R package jsonlite')
cfg <- jsonlite::fromJSON(args[3])
readtab <- function(p) read.delim(p,check.names=FALSE,stringsAsFactors=FALSE)
samples <- readtab(file.path(work,'samples.tsv'))
tpm <- readtab(file.path(work,'results/gene_TPM.tsv'))
rownames(tpm) <- tpm$gene_id;tpm$gene_id <- NULL;tpm <- as.matrix(tpm[,samples$sample,drop=FALSE])
metadata <- readtab(file.path(work,'results/gene_metadata.tsv'))
nr <- metadata$gene_id[tolower(as.character(metadata$is_nonreference)) %in% c('1','true')]
write_tsv <- function(x,p)write.table(x,p,sep='\t',row.names=FALSE,quote=FALSE)
plots <- isTRUE(cfg$plots)
figdir <- file.path(out,'figures');dir.create(figdir,showWarnings=FALSE)
if(plots){
  if(!requireNamespace('ggplot2',quietly=TRUE))stop('Install R package ggplot2')
  library(ggplot2)
  theme_set(theme_classic(base_size=11)+theme(plot.background=element_rect(fill='transparent',colour=NA),panel.background=element_rect(fill='transparent',colour=NA),legend.background=element_rect(fill='transparent',colour=NA),legend.key=element_rect(fill='transparent',colour=NA)))
  saveplot <- function(p,name,w=7,h=5){
    ggsave(file.path(figdir,paste0(name,'.png')),p,width=w,height=h,dpi=320,bg='transparent')
    ggsave(file.path(figdir,paste0(name,'.pdf')),p,width=w,height=h,device=cairo_pdf,bg='transparent')
  }
  labels <- if('group' %in% names(samples) && length(unique(samples$group))>1)samples$group else samples$tissue
  pca <- function(mat,name){
    x <- log2(mat+1);x <- x[apply(x,1,var)>0,,drop=FALSE]
    if(nrow(x)<2||ncol(x)<3){writeLines('Insufficient variable genes or samples for two PCA dimensions.',file.path(figdir,paste0(name,'.not_estimable.txt')));return(invisible(NULL))}
    fit <- prcomp(t(x),center=TRUE,scale.=FALSE);v <- fit$sdev^2/sum(fit$sdev^2)*100
    d <- data.frame(sample=samples$sample,group=labels,PC1=fit$x[,1],PC2=fit$x[,2])
    write_tsv(d,file.path(figdir,paste0(name,'.scores.tsv')))
    saveplot(ggplot(d,aes(PC1,PC2,color=group))+geom_point(size=3)+labs(x=sprintf('PC1 (%.1f%%)',v[1]),y=sprintf('PC2 (%.1f%%)',v[2]),color='Group',subtitle='Centered log2(TPM + 1); no gene variance scaling'),name)
  }
  pca(tpm,'PCA_all_genes_TPM');pca(tpm[intersect(nr,rownames(tpm)),,drop=FALSE],'PCA_nonreference_TPM')
  counts <- data.frame(sample=samples$sample,tissue=samples$tissue,genes=colSums(tpm[intersect(nr,rownames(tpm)),,drop=FALSE]>=1))
  write_tsv(counts,file.path(figdir,'nonreference_expression_per_sample.tsv'))
  saveplot(ggplot(counts,aes(tissue,genes,color=tissue))+geom_boxplot(outlier.shape=NA)+geom_point(position=position_jitter(width=.12,height=0,seed=42))+labs(x=NULL,y='Nonreference gene loci with TPM >= 1',color='Tissue'),'nonreference_expression_by_tissue')
  tissue <- readtab(file.path(out,'report/genes_tissue_mean_TPM_gt1.tsv'))
  for(which_set in c('all','nonreference')){
    d <- if(which_set=='all')tissue else tissue[tissue$is_nonreference==1,,drop=FALSE]
    if(!nrow(d))next
    sets <- split(d$gene_id,d$tissue);ids <- sort(unique(d$gene_id))
    membership <- vapply(sets,function(s)ids %in% s,logical(length(ids)))
    if(is.null(dim(membership)))membership <- matrix(membership,ncol=1,dimnames=list(NULL,names(sets)))
    pattern <- apply(membership,1,function(z)paste(as.integer(z),collapse=''))
    sizes <- sort(table(pattern),decreasing=TRUE)
    full <- data.frame(pattern=names(sizes),genes=as.integer(sizes));write_tsv(full,file.path(figdir,paste0('UpSet_',which_set,'.intersections.tsv')))
    top <- head(full,31);top$column <- seq_len(nrow(top))
    # Full exact intersections remain in TSV; dense views explicitly show top 31.
    title <- if(nrow(full)>31)'Top 31 exact intersections; full table exported' else 'All observed exact intersections'
    saveplot(ggplot(top,aes(factor(column),genes))+geom_col(fill='#4477AA')+labs(x='Exact intersection (see exported pattern table)',y='Unique gene loci',subtitle=title),paste0('UpSet_',which_set,'_intersection_sizes'),10,4)
  }
}
if(isTRUE(cfg$differential)){
  if(!requireNamespace('edgeR',quietly=TRUE)||!requireNamespace('statmod',quietly=TRUE))stop('Install R packages edgeR and statmod')
  control <- cfg$control;treatment <- cfg$treatment
  if(is.null(control)||is.null(treatment)||control==treatment)stop('Set distinct control and treatment groups')
  samples <- samples[samples$group %in% c(control,treatment),,drop=FALSE]
  if(any(table(factor(samples$group,levels=c(control,treatment)))<2))stop('At least two independent libraries per group required')
  groups <- factor(samples$group,levels=c(control,treatment));design <- model.matrix(~groups)
  if(nrow(design)<=ncol(design))stop('No residual degrees of freedom')
  raw <- readtab(file.path(work,'results/gene_counts.raw.tsv'));rownames(raw)<-raw$gene_id;raw$gene_id<-NULL
  raw <- as.matrix(raw[,samples$sample,drop=FALSE])
  if(any(!is.finite(raw))||any(raw<0)||any(raw!=round(raw)))stop('Differential analysis requires nonnegative integer raw counts')
  y <- edgeR::DGEList(raw,group=groups)
  y <- edgeR::calcNormFactors(y,method='TMM')
  minimum <- if(is.null(cfg$min_samples))2 else as.integer(cfg$min_samples)
  if(minimum<1||minimum>ncol(y))stop('Invalid analysis.min_samples')
  keep <- rowSums(edgeR::cpm(y)>=cfg$cpm)>=minimum
  dest <- file.path(out,'differential_expression');dir.create(dest,showWarnings=FALSE)
  write_tsv(data.frame(gene_id=rownames(y),retained=keep),file.path(dest,'expression_filter.tsv'))
  write_tsv(samples,file.path(dest,'analysis_samples.tsv'))
  if(sum(keep)<2)stop('Fewer than two genes pass count-expression filter; no differential model fitted')
  y <- y[keep,,keep.lib.sizes=FALSE];y <- edgeR::calcNormFactors(y,method='TMM')
  y <- edgeR::estimateDisp(y,design,robust=TRUE,trend.method=if(nrow(y)<50)'none' else 'locfit')
  fit <- edgeR::glmQLFit(y,design,robust=TRUE);test <- edgeR::glmQLFTest(fit,coef=2)
  tab <- edgeR::topTags(test,n=Inf,sort.by='PValue')$table
  tab <- data.frame(gene_id=rownames(tab),tab,is_nonreference=rownames(tab) %in% nr,row.names=NULL)
  tab$direction <- ifelse(tab$logFC>0,'Up','Down')
  tab$selected_FDR05_FC1 <- abs(tab$logFC)>1 & tab$FDR<.05
  tab$selected_P05_FC1 <- abs(tab$logFC)>1 & tab$PValue<.05
  write_tsv(tab,file.path(dest,'all_tested_genes.tsv'))
  for(rule in c('FDR05_FC1','P05_FC1')){
    selected <- tab[tab[[paste0('selected_',rule)]],,drop=FALSE]
    write_tsv(selected,file.path(dest,paste0(rule,'.all_selected.tsv')))
    selected <- selected[selected$is_nonreference,,drop=FALSE]
    write_tsv(selected,file.path(dest,paste0(rule,'.nonreference_selected.tsv')))
  }
  write_tsv(data.frame(sample=colnames(y),y$samples),file.path(dest,'library_normalization.tsv'))
  jsonlite::write_json(list(test='edgeR robust QL F-test',contrast=paste(treatment,'/',control),tested=nrow(tab),FDR05_FC1=sum(tab$selected_FDR05_FC1),P05_FC1=sum(tab$selected_P05_FC1),nonreference_FDR05_FC1=sum(tab$selected_FDR05_FC1 & tab$is_nonreference),nonreference_P05_FC1=sum(tab$selected_P05_FC1 & tab$is_nonreference)),file.path(dest,'summary.json'),auto_unbox=TRUE,pretty=TRUE)
  if(plots){
    for(rule in c('FDR05_FC1','P05_FC1')){
      tab$Expression <- ifelse(tab[[paste0('selected_',rule)]],tab$direction,'NS')
      tab$ordinate <- -log10(pmax(if(rule=='P05_FC1')tab$PValue else tab$FDR,.Machine$double.xmin))
      p <- ggplot(tab,aes(logFC,ordinate,color=Expression))+geom_point(size=1.1,alpha=.7)+scale_color_manual(values=c(Up='#36C989',Down='#9775CC',NS='#BBBBBB'))+geom_vline(xintercept=c(-1,1),linetype=2)+geom_hline(yintercept=-log10(.05),linetype=2)+geom_point(data=tab[tab$is_nonreference,,drop=FALSE],shape=5,color='black',size=2)+labs(x=paste0('log2FC (',treatment,' / ',control,')'),y=if(rule=='P05_FC1')'-log10(raw P)' else '-log10(FDR)',subtitle=if(rule=='P05_FC1')'Exploratory raw-P candidates; not FDR-controlled' else 'BH FDR-controlled selection')
      saveplot(p,paste0('Volcano_',rule))
    }
    selected <- tab$gene_id[tab$selected_P05_FC1 & tab$is_nonreference]
    if(length(selected)>=2){
      x <- log2(tpm[selected,samples$sample,drop=FALSE]+1)
      fit <- cmdscale(dist(t(x)),k=2,eig=TRUE)
      scores <- data.frame(sample=samples$sample,group=groups,Dimension1=fit$points[,1],Dimension2=fit$points[,2])
      write_tsv(scores,file.path(figdir,'MDS_nonreference_rawP_candidates_TPM.coordinates.tsv'))
      saveplot(ggplot(scores,aes(Dimension1,Dimension2,color=group))+geom_point(size=3)+labs(color='Stage',subtitle='Descriptive MDS after candidate selection; log2(TPM + 1)'),'MDS_nonreference_rawP_candidates_TPM')
    }else writeLines('Fewer than two selected nonreference candidates; MDS not estimated.',file.path(figdir,'MDS_not_estimable.txt'))
  }
}
writeLines(capture.output(sessionInfo()),file.path(out,'report/R_sessionInfo.txt'))

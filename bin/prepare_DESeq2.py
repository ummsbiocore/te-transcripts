#!/usr/bin/env python3

from collections import defaultdict
from textwrap import dedent
import os
import sys
import argparse

def main(args):

	comparisons = read_comparisons(args.comparisons)
	highlighted = ', '.join(["'%s'" % i for i in args.highlighted_genes])

	for comparison in comparisons:
		process_comparison(comparison[0], comparison[1], comparison[2], comparison[3], highlighted, args)

def read_comparisons(comparison_infile):

	comparisons = []

	header = open(comparison_infile).readline().rstrip('\n').split('\t')
	if len(header) > 1:
		sep = '\t'
	else:
		sep = ','

	with open(comparison_infile) as infile:
		infile.readline()
		for line in infile:
			cur = line.rstrip('\n').split(sep)

			if len(cur) == 3:
				comparisons.append((cur[1], cur[0], cur[2], 'group'))
			else:
				if cur[3] == '':
					comparisons.append((cur[1], cur[0], cur[2], 'group'))
				else:
					comparisons.append((cur[1], cur[0], cur[2], cur[3]))
	
	return comparisons

def read_groups(metadata_infile, column):

	groups = defaultdict(list)

	header = open(metadata_infile).readline().rstrip('\n').split('\t')
	if len(header) > 1:
		sep = '\t'
	else:
		sep = ','

	header = open(metadata_infile).readline().rstrip('\n').split(sep)

	if column in header:
		index = header.index(column)
	else:
		index = 1

	with open(metadata_infile) as infile:
		infile.readline()
		for line in infile:
			cur = line.rstrip('\n').split(sep)
			groups[cur[index]].append(cur[0])

	return groups

def process_comparison(treatment, control, name, column, highlighted, args):

	groups = read_groups(args.groups, column)

	if args.design.replace('~', '') != column:
		design = '~%s' % column
	else:
		design = args.design

	build_report(groups, design, treatment, control, name, column, highlighted, args, '')
	
	if args.excluded_events != '':
		build_report(groups, design, treatment, control, '%s_excludelist' % (name), column, highlighted, args, args.excluded_events)

def build_report(groups, design, treatment, control, name, column, highlighted, args, excluded_events):
	output = []
	output.append(print_header(name))
	output.append(print_libraries())
	output.append(print_functions())
	output.append(choose_modules(args))
	output.append(read_counts(groups, treatment, control, args, excluded_events))
	output.append(filter_counts(args))
	output.append(batch_correction(args))
	output.append(sample_info())
	output.append(quality_control())
	output.append(count_distribution(column))
	output.append(all2all())
	output.append(run_pca(args))
	output.append(run_deseq(treatment, control, highlighted, args, name, design, column))
	output.append(volcano_plot())
	output.append(ma_plot())
	output.append(heatmap())
	output.append(write_tables(name))
	output.append(write_warnings())
	output.append(session_info())
	write_output(name, output)

def print_header(name):

	return(dedent(
	'''
	---
	title: {}
	date: "`r Sys.Date()`"
	output:
	  html_document:
	    code_folding: hide
	---
	'''.format(name)).strip())

def print_libraries():
	
	return(dedent(
	'''
	```{r, load_libraries, message=FALSE, include=FALSE}
	# Load Libraries
	library(dplyr)
	library(ggplot2)
	library(tidyr)
	library(DESeq2)
	library(scales)
	library(ggrepel)
	library(clusterProfiler)
	library(sva)
	library(gplots)
	library(yaml)
	library(DT)
	library(htmltools)
	```'''))

def print_functions():

	return(dedent(
	'''
	```{r local_functions, include=FALSE}
	# Load custom functions
	reverselog = function() {
	  trans_new("reverselog", function(x) -log10(x), function(x) 10^(-x), log_breaks(base = 10), domain = c(1e-1000, Inf))
	}
	
	quiet <- function(x) { 
	  sink(tempfile()) 
	  on.exit(sink()) 
	  invisible(force(x)) 
	} 
	
	read_metadata = function(metadata_file) {
	  metadata = read.delim(metadata_file, header=TRUE)
	  if (length(metadata) < 2) {
	    metadata = read.delim(metadata_file, sep=',', header=TRUE)
	  }
	  return(metadata)
	}
	
	read_counts = function(count_file, samples, excluded_events) {
	  
	  counts = read.delim(count_file, check.names=FALSE)
	  if (length(counts) < 2) {
	    counts = read.delim(count_file, check.names=FALSE, sep=',')
	  }

	  counts = counts %>%
	           select(gene, all_of(samples)) %>%
	           filter(!(gene %in% excluded_events)) %>%
	           tibble::column_to_rownames('gene')
	  counts = as.matrix(counts)
	  mode(counts) = 'integer'
	  return(counts)
	}
	
	check_samples = function(all_counts, total_count_min) {
	  bad_samples = (data.frame(Counts=colSums(all_counts)) %>% tibble::rownames_to_column('Sample') %>% filter(Counts < total_count_min))$Sample
	  return(bad_samples)
	}
	
	sample_filter = function(all_counts, bad_samples) {
	  x = as.data.frame(all_counts) %>% select(-all_of(bad_samples))
	  return(as.matrix(x))
	}

	count_distribution = function(counts, samples, min_counts_per_event, group_by) {
	
	  options(scipen=999)
	
	  needed_colors = nrow(samples %>% distinct(!!sym(group_by)))
	  colors = c("#4682b4", "#000000", "#6d2578", "#DC6D00", "#1e6234", "#FFFF33", "#A65628", "#F781BF", "#999999")

	  df = as.data.frame(counts) %>% tibble::rownames_to_column("gene") %>%
	       pivot_longer(!gene, names_to='sample_name', values_to='count') %>%
	       left_join(samples, by='sample_name') %>%
	       group_by(gene, !!sym(group_by)) %>%
	       summarise(ave_count = mean(count), .groups='drop') %>%
	       filter(ave_count > 0)
	
	  ggplot(df, aes(x=ave_count, fill=as.factor(!!sym(group_by)))) +
	    theme_classic() +
	    theme(strip.background = element_blank(),
	          strip.text = element_blank()) +
	    facet_wrap(as.formula(paste("~", group_by)), ncol=1) +
	    scale_x_continuous(trans='log10', breaks=c(.1,1,10,100,1000,10000,100000,1000000), name='Raw Counts') +
	    scale_y_continuous(name='Number of Genes', expand=c(0,0)) +
	    {if (needed_colors <= length(colors)) scale_fill_manual(values = colors)} +
	    geom_histogram(bins=100) +
	    geom_vline(xintercept = min_counts_per_event, linetype=2, color='firebrick')
	}
	
	filter_counts = function(input, min_counts_per_event, min_samples_per_event) {
	  keep = rowSums(input >= min_counts_per_event) >= min_samples_per_event
	  return(input[keep,])
	}
	
	run_pca = function(input, retx=TRUE, center=TRUE, scale=TRUE, transformation='Default') {
	  if (transformation == 'None') {
	    keep = subset(input, apply(input, 1, var, na.rm = TRUE) >  0)
	    return(prcomp(t(keep), retx = retx, center = center, scale. = scale))
	  } else if (transformation == 'vst' | (transformation == 'Default' & ncol(input) > 50)) {
	    transformed = vst(input)
	    return(prcomp(t(transformed), retx = retx, center = center, scale. = FALSE))
	  } else {
	    transformed = rlog(input)
	    return(prcomp(t(transformed), retx = retx, center = center, scale. = FALSE))
	  }
	}
	
	batch_correction = function(df, batch_norm_method, batch_correction_column, batch_correction_group_column) {
	  set.seed(1)
	  norm_df = getNormalizedMatrix(df, method=batch_norm_method)
	  mode(norm_df) = 'integer'
	  adjusted_df = apply(norm_df, 2, function(x) as.integer(x) + runif(1, 0, 0.01))
	  colnames(adjusted_df) <- colnames(norm_df)
	  rownames(adjusted_df) <- rownames(norm_df)
	  corrected = quiet(ComBat_seq(adjusted_df, batch = batch_correction_column, group=batch_correction_group_column))
	  mode(corrected) = 'integer'
	  return(corrected)
	}	

	call_significance = function(df, padj_significance_cutoff, fc_significance_cutoff, padj_floor, fc_ceiling, num_labeled, apply_shrinkage) {
	  labels = df %>%
	           filter(padj < padj_significance_cutoff & abs(log2FoldChange) > fc_significance_cutoff) %>%
	           arrange(-abs(log2FoldChange)) %>%
	           mutate(Rank = row_number()) %>%
	           filter(Rank <= num_labeled) %>%
	           select(gene)
	
	  if (apply_shrinkage) {
	    df = df %>% 
	         mutate(Direction = case_when(padj < padj_significance_cutoff & log2FoldChange_shrink > fc_significance_cutoff ~ 'Upregulated',
	                                      padj < padj_significance_cutoff & log2FoldChange_shrink < -fc_significance_cutoff ~ 'Downregulated',
	                                      TRUE ~ 'No Change')) %>%
	         mutate(log2FoldChange_shrink = case_when(log2FoldChange_shrink < 0 & abs(log2FoldChange_shrink) > fc_ceiling ~ -fc_ceiling,
	                                                  log2FoldChange_shrink > 0 & abs(log2FoldChange_shrink) > fc_ceiling ~ fc_ceiling,
	                                                  TRUE ~ log2FoldChange_shrink))
	  } else {
	    df = df %>% 
	         mutate(Direction = case_when(padj < padj_significance_cutoff & log2FoldChange > fc_significance_cutoff ~ 'Upregulated',
	                                      padj < padj_significance_cutoff & log2FoldChange < -fc_significance_cutoff ~ 'Downregulated',
	                                      TRUE ~ 'No Change'))
	  }
	
	  df = df %>%
	       mutate(Direction = factor(Direction, levels=c("Upregulated", "No Change", "Downregulated"))) %>%
	       mutate(Significant = case_when(Direction == 'No Change' ~ 'Non-Significant',
	                                      TRUE ~ 'Significant')) %>%
	       mutate(padj = case_when(padj < padj_floor ~ padj_floor,
	                               TRUE ~ padj)) %>%
	       mutate(log2FoldChange = case_when(log2FoldChange < 0 & abs(log2FoldChange) > fc_ceiling ~ -fc_ceiling,
	                                         log2FoldChange > 0 & abs(log2FoldChange) > fc_ceiling ~ fc_ceiling,
	                                         TRUE ~ log2FoldChange)
	       ) %>%
	       mutate(Label = case_when(((Significant == 'Significant' & gene %in% labels$gene) | Highlighted == TRUE) ~ alias)) %>%
	       mutate(Group = case_when(Highlighted==TRUE & Direction=='Upregulated' ~ 'Upregulated_Highlighted',
	                                Highlighted==TRUE & Direction=='No Change' ~ 'NoChange_Highlighted',
	                                Highlighted==TRUE & Direction=='Downregulated' ~ 'Downregulated_Highlighted',
	                                Highlighted==FALSE & Direction=='Upregulated' ~ 'Upregulated',
	                                Highlighted==FALSE & Direction=='No Change' ~ 'No Change',
	                                Highlighted==FALSE & Direction=='Downregulated' ~ 'Downregulated'
	                      )
	       )
	  return(df)
	}
	
	add_alias = function(df, add_alias=TRUE, fromType='ENSEMBL', toType='SYMBOL', org='org.Hs.eg.db') {
	  if (add_alias) {
	    BiocManager::install('org.Hs.eg.db', update=FALSE, force=TRUE)
	    mapped_names <- data.frame(bitr(df$gene, fromType=fromType, toType=toType, OrgDb='org.Hs.eg.db'))
	  
	    uniquely_mapped = mapped_names %>% 
	                      group_by(!!sym(toType)) %>% 
	                      summarise(count = n()) %>% 
	                      filter(count == 1) %>%
	                     select(-count)
	    colnames(uniquely_mapped) = c('alias')
	  
	    mapped_names = mapped_names %>% filter(!!sym(toType) %in% uniquely_mapped$alias)
	  
	    return(df %>%
	           left_join(mapped_names, by=c(gene = fromType)) %>%
	           mutate(alias = case_when(!is.na(!!sym(toType)) ~ !!sym(toType),
	                                    TRUE ~ gene)) %>%
	           select(-!!sym(toType))
	    )
	    } else {
	      return(df %>% mutate(alias=gene))
	    }
	}
	
	add_highlights = function(df, to_highlight) {
	  return(
	    df %>% 
	    mutate(Highlighted = 
	      case_when((alias %in% to_highlight | gene %in% to_highlight) ~ TRUE,
	                 TRUE ~ FALSE
	      )
	    )
	  )
	}
	
	post_processing = function(deseq_res, padj_significance_cutoff, fc_significance_cutoff, num_labeled, highlighted, add_alias=TRUE, fromType='ENSEMBL', toType='SYMBOL', org='org.Hs.eg.db', apply_shrinkage=FALSE) {
	  post = add_alias(deseq_res, add_alias=add_alias, fromType=fromType, toType=toType, org='org.Hs.eg.db')
	  post = add_highlights(post, highlighted)
	  post = call_significance(post, padj_significance_cutoff, fc_significance_cutoff, padj_floor, fc_ceiling, num_labeled, apply_shrinkage)
	  return(post)
	}
	
	pca_plot = function(pca, metadata, x='PC1', y='PC2', color_by='', shape_by='', fill_by='', alpha_by='', label_by='', size=4) {
	
	  col_names <- names(metadata)
	  metadata[,col_names] <- lapply(metadata[,col_names] , factor)
	
	  if (color_by=='NA') { color_by='' }
	  if (shape_by=='NA') { shape_by='' }
	  if (fill_by=='NA') { fill_by='' }
	  if (alpha_by=='NA') { alpha_by='' }
	  if (label_by=='NA') { label_by='' }
	  
	  if (color_by != '') {
	    needed_colors = nrow(metadata %>% distinct(!!sym(color_by)))
	  } else {
	    needed_colors = Inf
	  }
	  colors = c("#b22222", "#4682b4", "#6d2578", "#DC6D00", "#1e6234", "#FFFF33", "#A65628", "#F781BF", "#999999")

	  
	  PoV = data.frame(PoV=pca$sdev^2/sum(pca$sdev^2)*100) %>% mutate(PC = paste0("PC", row_number()))
	  pca_results = as.data.frame(pca$x) %>%
	                tibble::rownames_to_column('sample_name') %>%
	                left_join(metadata, by='sample_name')
	  
	  PoV_X = round((PoV %>% filter(PC==x))$PoV, 2)
	  PoV_Y = round((PoV %>% filter(PC==y))$PoV, 2)
	
	  ggplot(pca_results, aes(x=!!sym(x), y=!!sym(y), color=!!sym(color_by), shape=!!sym(shape_by), fill=!!sym(fill_by), alpha=!!sym(alpha_by), label=!!sym(label_by))) +
	    theme_classic() +
	    xlab(paste0(x, ': ', PoV_X, '%')) +
	    ylab(paste0(y, ': ', PoV_Y, '%')) +
	    {if (needed_colors <= length(colors)) scale_color_manual(values = colors)} +
	    geom_point(size=size) +
	    {if (label_by != '') geom_label(label.size=NA, fill=NA, color='black')}
	}
	
	scree_plot = function(pca) {
	
	PoV = data.frame(PoV=pca$sdev^2/sum(pca$sdev^2)*100) %>%
	      mutate(PC = row_number()) %>%
	      mutate(Label = case_when(PoV >=.01 ~ as.character(round(PoV, 2)),
	                               TRUE ~ "<.01"))
	
	ggplot(PoV, aes(x=PC, y=PoV, label=Label)) +
	  theme_classic() +
	  scale_x_continuous(breaks=seq(1, min(nrow(PoV), 20), 1), limits=c(0, min(nrow(PoV), 20))) +
	  scale_y_continuous(name="Percent of Variation", limits=c(0,100), expand=c(0,0)) +
	  geom_bar(stat='identity', fill='black') +
	  geom_label(fill=NA, label.size = NA, vjust=-.05)
	}
	
	ma_plot = function(deseq_res, Y='log2FoldChange', padj_cutoff=.05, fc_cutoff=1, padj_floor=0, fc_ceiling=Inf, num_labeled=Inf,
	                  upregulated_color='firebrick', noChange_color='grey', downregulated_color='steelblue',
	                  upregulated_highlight_color='#1e6234', noChange_highlight_color='black', downregulated_highlight_color='#1e6234',
	                  upregulated_alpha=1, noChange_alpha=.3, downregulated_alpha=1,
	                  upregulated_highlight_alpha=1, noChange_highlight_alpha=1, downregulated_highlight_alpha=1,
	                  upregulated_size=2, noChange_size=1, downregulated_size=2,
	                  upregulated_highlight_size=2, noChange_highlight_size=2, downregulated_highlight_size=2,
	                  fc_markers=TRUE, fc_marker_color='grey', center_marker=TRUE, center_marker_color='black',
	                  display_upregulated_count = TRUE, display_downregulated_count = TRUE, display_noChange_count=FALSE) {
	
	  options(scipen=999)
	
	  colors = c('Upregulated'=upregulated_color, 'No Change'=noChange_color, 'Downregulated'=downregulated_color, 'Upregulated_Highlighted'=upregulated_highlight_color, 'NoChange_Highlighted'=noChange_highlight_color, 'Downregulated_Highlighted'=downregulated_highlight_color)
	  alphas = c('Upregulated'=upregulated_alpha, 'No Change'=noChange_alpha, 'Downregulated'=downregulated_alpha, 'Upregulated_Highlighted'=upregulated_highlight_alpha, 'NoChange_Highlighted'=noChange_highlight_alpha, 'Downregulated_Highlighted'=downregulated_highlight_alpha)
	  sizes = c('Upregulated'=upregulated_size, 'No Change'=noChange_size, 'Downregulated'=downregulated_size, 'Upregulated_Highlighted'=upregulated_highlight_size, 'NoChange_Highlighted'=noChange_highlight_size, 'Downregulated_Highlighted'=downregulated_highlight_size)
	
	  upregulated_count = nrow(deseq_res %>% filter(Direction == 'Upregulated'))
	  noChange_count = nrow(deseq_res %>% filter(Direction == 'No Change'))
	  downregulated_count = nrow(deseq_res %>% filter(Direction == 'Downregulated'))
	
	  return(
	    ggplot(deseq_res, aes(x=baseMean, y=!!sym(Y), color=Group, alpha=Group, size=Group, label=Label)) +
	      theme_classic() +
	      theme(legend.position = 'none') +
	      scale_x_continuous(trans='log10', name="Expression") +
	      scale_y_continuous(name="Fold Change (log2)") +
	      scale_color_manual(values=colors) +
	      scale_alpha_manual(values=alphas) +
	      scale_size_manual(values=sizes) +
	      geom_point() +
	      {if (fc_markers) geom_hline(yintercept = fc_cutoff, linetype=2, color=fc_marker_color)} +
	      {if (center_marker) geom_hline(yintercept = 0, linetype=2, color=center_marker_color)} +
	      {if (fc_markers) geom_hline(yintercept = -fc_cutoff, linetype=2, color=fc_marker_color)} +
	      {if (display_upregulated_count) annotate('text', x=Inf, y=Inf, hjust=1, vjust=1, color=upregulated_color, label=paste0("Upregulated: ", upregulated_count))} +
	      {if (display_downregulated_count) annotate('text', x=Inf, y=-Inf, hjust=1, vjust=-1, color=downregulated_color, label=paste0("Downregulated: ", downregulated_count))} +
	      {if (display_noChange_count) annotate('text', x=Inf, y=fc_cutoff, hjust=1, vjust=1.2, color=noChange_color, label=paste0("No Change: ", noChange_count))} +
	      geom_label_repel(label.size=NA, fill=NA, na.rm=TRUE, max.overlaps = 50, max.time = 5)
	  )
	}
	
	volcano_plot = function(deseq_res, X='log2FoldChange', padj_cutoff=.05, fc_cutoff=1, padj_floor=0, fc_ceiling=Inf, num_labeled=Inf,
	                        upregulated_color='firebrick', noChange_color='grey', downregulated_color='steelblue',
	                        upregulated_highlight_color='#1e6234', noChange_highlight_color='black', downregulated_highlight_color='#1e6234',
	                        upregulated_alpha=1, noChange_alpha=.3, downregulated_alpha=1,
	                        upregulated_highlight_alpha=1, noChange_highlight_alpha=1, downregulated_highlight_alpha=1,
	                        upregulated_size=2, noChange_size=1, downregulated_size=2,
	                        upregulated_highlight_size=2, noChange_highlight_size=2, downregulated_highlight_size=2,
	                        fc_markers=TRUE, fc_marker_color='grey', center_marker=TRUE, center_marker_color='black',
	                        padj_marker=TRUE, padj_marker_color='grey',
	                        display_upregulated_count = TRUE, display_downregulated_count = TRUE, display_noChange_count=FALSE) {
	
	  colors = c('Upregulated'=upregulated_color, 'No Change'=noChange_color, 'Downregulated'=downregulated_color, 'Upregulated_Highlighted'=upregulated_highlight_color, 'NoChange_Highlighted'=noChange_highlight_color, 'Downregulated_Highlighted'=downregulated_highlight_color)
	  alphas = c('Upregulated'=upregulated_alpha, 'No Change'=noChange_alpha, 'Downregulated'=downregulated_alpha, 'Upregulated_Highlighted'=upregulated_highlight_alpha, 'NoChange_Highlighted'=noChange_highlight_alpha, 'Downregulated_Highlighted'=downregulated_highlight_alpha)
	  sizes = c('Upregulated'=upregulated_size, 'No Change'=noChange_size, 'Downregulated'=downregulated_size, 'Upregulated_Highlighted'=upregulated_highlight_size, 'NoChange_Highlighted'=noChange_highlight_size, 'Downregulated_Highlighted'=downregulated_highlight_size)
	
	  upregulated_count = nrow(deseq_res %>% filter(Direction == 'Upregulated'))
	  noChange_count = nrow(deseq_res %>% filter(Direction == 'No Change'))
	  downregulated_count = nrow(deseq_res %>% filter(Direction == 'Downregulated'))
	
	  return(
	    ggplot(deseq_res, aes(x=!!sym(X), y=padj, color=Group, alpha=Group, size=Group, label=Label)) +
	      theme_classic() +
	      theme(legend.position = 'none') +
	      scale_x_continuous(name="Fold Change (log2)") +
	      scale_y_continuous(trans=reverselog(), name='Significance', labels=trans_format('log10',math_format(10^.x))) +
	      scale_color_manual(values=colors) +
	      scale_alpha_manual(values=alphas) +
	      scale_size_manual(values=sizes) +
	      geom_point() +
	      {if (fc_markers) geom_vline(xintercept = fc_cutoff, linetype=2, color=fc_marker_color)} +
	      {if (center_marker) geom_vline(xintercept = 0, linetype=2, color=center_marker_color)} +
	      {if (fc_markers) geom_vline(xintercept = -fc_cutoff, linetype=2, color=fc_marker_color)} +
	      {if (padj_marker) geom_hline(yintercept = padj_cutoff, linetype=2, color=padj_marker_color)} +
	      {if (display_upregulated_count) annotate('text', x=Inf, y=0, hjust=1, vjust=1, color=upregulated_color, label=paste0("Upregulated: ", upregulated_count))} +
	      {if (display_downregulated_count) annotate('text', x=-Inf, y=0, hjust=-.01, vjust=1, color=downregulated_color, label=paste0("Downregulated: ", downregulated_count))} +
	      {if (display_noChange_count) annotate('text', x=0, y=0, vjust=1, color=noChange_color, label=paste0("No Change: ", noChange_count))} +
	      geom_label_repel(label.size=NA, fill=NA, na.rm=TRUE, max.overlaps = 50, max.time = 5)
	  )
	}
	
	heatmap_plot = function(df) {
	  norm = getNormalizedMatrix(df, method='TMM')
	  ld = log2(norm+.1)
	  cldt = scale(t(ld), center=TRUE, scale=TRUE)
	  cld = t(cldt)
	  
	  if (nrow(cld) > 1) {
	  	heatmap.2(cld, Rowv=TRUE, dendrogram='column', Colv=TRUE, col=bluered(256), 
	              labRow=NA, density.info='none', trace='none', cexCol=.8,
	              hclust=function(x) hclust(x,method="complete"),
	              distfun=function(x) as.dist((1-cor(t(x)))/2)
	    )
	  } else {
	      x = as.data.frame(cld) %>% tibble::rownames_to_column('Gene') %>%
	          pivot_longer(!Gene, names_to = "Sample", values_to = "Value")
	
	      ggplot(x, aes(x=Sample, y=Gene, fill=Value)) +
	             theme_classic() +
	             theme(axis.line = element_blank(),
	                   axis.ticks= element_blank()) +
	             scale_fill_gradient(low="#0000FF" , high='#FF0000') +
	             geom_tile()
	  }
	}
	
	all2all <- function(data, cex=2) {
	    pcor <- function(x, y, ...) panel.cor(x, y, cex.cor = cex)
	    nr <- nrow(data)
	    if (nr > 1000)
	        nr <- 1000
	    pairs(log10(data[1:nr, ]), cex = 0.25,
	            diag.panel = panel.hist, lower.panel = pcor)
	}

	panel.hist <- function(x, ...) {
	    usr <- par("usr")
	    on.exit(par(usr))
	    par(usr = c(usr[1:2], 0, 1.5))
	    h <- hist(x, plot = FALSE)
	    breaks <- h$breaks
	    nb <- length(breaks)
	    y <- h$counts
	    y <- y / max(y)
	    rect(breaks[-nb], 0, breaks[-1], y, col = "red", ...)
	}
	
	panel.cor <- function(x, y, prefix = "rho=", cex.cor=2, ...){
	    usr <- par("usr")
	    on.exit(par(usr))
	    par(usr = c(0, 1, 0, 1))
	    r <- cor.test(x, y, method = "spearman",
	        na.rm = TRUE, exact = FALSE)$estimate
	    txt <- round(r, digits = 2)
	    txt <- paste0(prefix, txt)
	    text(0.5, 0.5, txt, cex = cex.cor)
	}

	getNormalizedMatrix <- function(M = NULL, method = "TMM") {
	    if (is.null(M) ) return (NULL)
	    M[is.na(M)] <- 0
	    norm <- M
	    if (!(method == "none" || method == "MRN")){
	        norm.factors <- edgeR::calcNormFactors(M, method = method)
	        norm <- edgeR::equalizeLibSizes(edgeR::DGEList(M,
	            norm.factors = norm.factors))$pseudo.counts
	    }else if(method == "MRN"){
	        columns <- colnames(M)
	        conds <- columns
	        coldata <- prepGroup(conds, columns)
	        M[, columns] <- apply(M[, columns], 2,
	            function(x) as.integer(x))
	        dds <- DESeqDataSetFromMatrix(countData = as.matrix(M),
	            colData = coldata, design = ~group)
	        dds <- estimateSizeFactors(dds)
	        norm <- counts(dds, normalized=TRUE)
	    }
	    return(norm)
	}

	deseq_results_table = function(post_res, apply_shrinkage, name) {	
	
	  if (apply_shrinkage) {
	  
	    sketch = htmltools::withTags(table(
	      class = 'display',
	      thead(
	        tr(
	          th(rowspan=2, 'gene'), th(rowspan=2, 'baseMean'), th(rowspan=2, 'pvalue'), th(rowspan=2, 'padj'), th(colspan=3, style = "text-align:center; border-left: solid 1px;", "Standard"), th(colspan=2, style = "text-align:center; border-left: solid 1px;", "Shrinked"), th(rowspan=2, 'Direction'), th(rowspan=2, 'Rank'), th(rowspan=2, '-|log2FoldChange|'),
	        ),
	        tr(
	          th("log2FoldChange", style="border-left: solid 1px;"), th("lfcSE"), th("stat"), th("log2FoldChange", style="border-left: solid 1px;"), th("lfcSE")
	        )
	      )
	    ))
	
	    rank = post_res %>% filter(Significant=='Significant') %>% mutate(Rank = row_number(abs(log2FoldChange_shrink))) %>% select(gene, Rank)
	    df = post_res %>% left_join(rank, by='gene') %>% select(gene, baseMean, pvalue, padj, log2FoldChange, lfcSE, stat, log2FoldChange_shrink, lfcSE_shrink, Direction, Rank) %>% mutate(abs_fc=-abs(log2FoldChange_shrink)) %>% arrange(Rank, abs_fc)
	    return(
	      datatable(df,
	        rownames=FALSE,
	        extensions = 'Buttons',
	        options=list(
	          columnDefs=list(list(visible=FALSE, targets=c(9, 10, 11))),
	          dom = 'lftBipr',
	          buttons = list(
	            list(extend = 'csvHtml5', text='Download', filename = paste0(name, "_deseq2_results"), extension='.tsv', fieldBoundary='', fieldSeparator='\t')
	           )
	          ),
	        container=sketch
	      ) %>% 
	      formatRound(columns=c('baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj', 'log2FoldChange_shrink', 'lfcSE_shrink'), digits=4) %>%
	      formatSignif(columns=c("pvalue", "padj"), digits=4) %>%
	      formatStyle(c(5, 8), `border-left` = "solid 1px") %>%
	      formatStyle('Direction',target = 'row', color = styleEqual(c("No Change", "Upregulated", "Downregulated"), c('black', 'firebrick', 'steelblue')))
	    )
	  } else {
	
	    rank = post_res %>% filter(Significant=='Significant') %>% mutate(Rank = row_number(abs(log2FoldChange))) %>% select(gene, Rank)
	    df = post_res %>% left_join(rank, by='gene') %>% select(gene, baseMean, pvalue, padj, log2FoldChange, lfcSE, stat, Direction, Rank) %>% mutate(abs_fc=-abs(log2FoldChange)) %>% arrange(Rank, abs_fc)
	
	    return(
	      datatable(df,
	        rownames=FALSE,
	        extensions = 'Buttons',
	        options=list(
	          columnDefs=list(list(visible=FALSE, targets=c(7,8,9))),
	          dom = 'lftBipr',
	          buttons = list(
	            list(extend = 'csvHtml5', text='Download', filename = paste0(name, "_deseq2_results"), extension='.tsv', fieldBoundary='', fieldSeparator='\t')
	           )
	          )
	      ) %>% 
	      formatRound(columns=c('baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj'), digits=4) %>%
	      formatSignif(columns=c("pvalue", "padj"), digits=4) %>%
	      formatStyle('Direction',target = 'row', color = styleEqual(c("No Change", "Upregulated", "Downregulated"), c('black', 'firebrick', 'steelblue')))
	    )
	  }
	}

	output_parameters = function(prefix, min_counts_per_event, min_samples_per_event, transformation, include_batch_correction, batch_correction_column, batch_correction_group_column, batch_norm_method, design, fitType, apply_shrinkage, shrinkage_type, padj_significance_cutoff, fc_significance_cutoff, padj_floor, fc_ceiling) {
	  params = data.frame(min_counts_per_event = min_counts_per_event,
	                      min_samples_per_event = min_samples_per_event,
	                      transformation = transformation,
	                      include_batch_correction = as.character(include_batch_correction),
	                      batch_correction_column = batch_correction_column,
	                      batch_correction_group_column = batch_correction_group_column,
	                      batch_norm_method = batch_norm_method,
	                      design=format(design),
	                      fitType=fitType,
	                      apply_shrinkage=apply_shrinkage,
	                      shrinkage_type=shrinkage_type,
	                      padj_significance_cutoff=padj_significance_cutoff,
	                      fc_significance_cutoff=fc_significance_cutoff,
	                      padj_floor=padj_floor,
	                      fc_ceiling=fc_ceiling
	           )
	  write_yaml(params, file=paste0('outputs/', prefix, '_params.yaml'))
	}
	```
	'''
	))

def choose_modules(args):

	return(dedent(
	'''
	```{{r, choose_modules}}
	include_count_distribution = {}
	include_all2all = {}
	include_pca = {}
	include_batch_correction = {}
	include_deseq2 = {}
	apply_shrinkage = {}
	convert_names = {}
	include_ma = {}
	include_volcano = {} 
	include_heatmap = {}
	include_volcano_highlighted = {}
	include_ma_highlighted = {}
	```'''.format(args.include_distribution, args.include_all2all, args.include_pca, args.include_batch_correction, args.include_DESeq2, args.apply_shrinkage, args.convert_names, args.include_ma, args.include_volcano, args.include_heatmap, args.include_volcano_highlighted, args.include_ma_highlighted)))

def read_counts(groups, treatment, control, args, excluded_events):

	treatments = ', '.join(["'%s'" % sample for sample in groups[treatment]])
	controls = ', '.join(["'%s'" % sample for sample in groups[control]])

	return(dedent(
	'''
	```{{r, load_counts}}
	excluded_events = c({})

	# Load count data
	all_metadata = read_metadata("{}{}")
	all_counts = read_counts("{}{}", all_metadata$sample_name, excluded_events)

	min_counts_per_sample = {}
	min_counts_per_event = {}
	min_samples_per_event = {}

	bad_samples = check_samples(all_counts, min_counts_per_sample)
	all_counts = sample_filter(all_counts, bad_samples)

	treatments = c({})
	good_treatments = setdiff(treatments, bad_samples)

	controls = c({})
	good_controls = setdiff(controls, bad_samples)
	
	selected_counts = read_counts("{}{}", c(good_treatments, good_controls), excluded_events)
	selected_metadata = data.frame(sample_name = colnames(selected_counts)) %>% left_join(all_metadata, by='sample_name')
	```'''.format(','.join('"%s"' % (i) for i in excluded_events), args.input_prefix, args.groups, args.input_prefix, args.counts, args.min_counts_per_sample, args.min_counts_per_event, args.min_samples_per_event, treatments, controls, args.input_prefix, args.counts)))

def filter_counts(args):

	return(dedent(
	'''
	```{{r, filter, echo=FALSE}}
	filter_type = '{}'

	all_counts_filtered = filter_counts(all_counts, min_counts_per_event, min_samples_per_event)
	
	if (filter_type=='local') {{
	  selected_counts_filtered = filter_counts(selected_counts, min_counts_per_event, min_samples_per_event) %>% as.matrix()
	}} else {{
	  selected_counts_filtered = as.data.frame(all_counts_filtered) %>% select(all_of(good_treatments), all_of(good_controls)) %>% as.matrix()
	}}
	```'''.format(args.filter_type)))

def batch_correction(args):

	return(dedent(
	'''
	```{{r, batch_correction, eval=include_batch_correction, echo=FALSE}}
	batch_norm_method = '{}'
	batch_correction_column = as.factor(all_metadata${})
	batch_correction_group_column = as.factor(all_metadata${})
	
	all_counts_batch_corrected = batch_correction(all_counts_filtered, batch_norm_method, batch_correction_column, batch_correction_group_column)
	
	selected_counts_batch_corrected = as.data.frame(all_counts_batch_corrected) %>% 
	                                  select(all_of(good_treatments), all_of(good_controls)) %>% 
	                                  tibble::rownames_to_column("Event") %>%
	                                  filter(Event %in% rownames(selected_counts_filtered)) %>%
	                                  tibble::column_to_rownames("Event") %>%
	                                  as.matrix()
	```'''.format(args.batch_correction_normalization_algorithm, args.batch_correction_column, args.batch_correction_group_column)))


def sample_info():

	return(dedent(
	'''
	```{r, sample_info, echo=FALSE, results='asis'}
	cat("## Sample Info")
	datatable(selected_metadata,
	  rownames = FALSE					
	)
	```'''))

def quality_control():

	return(dedent(
	'''
	```{r, QC, eval=include_count_distribution | include_all2all | include_pca, echo=FALSE, results='asis'}
	cat('## Quality Control {.tabset .tabset-pills}')
	```'''))

def count_distribution(column):

	return(dedent(
	'''
	```{{r, count_distribution, eval=include_count_distribution, echo=FALSE, results='asis'}}
	cat('### Raw Count Distribution\\n')
	cat('In Eukaryotes only a subset of all genes are expressed in a given cell. Expression is therefore a bimodal distribution, with non-expressed genes having counts that result from experimental and biological noise. It is important to filter out the genes that are not expressed before doing differential gene expression. You can decide which cutoff separates expressed vs non-expressed genes by looking at the following histogram. The current count threshold is at least ', min_counts_per_event, ' raw reads (red line) in at least ', min_samples_per_event, ' samples.\\n')
	
	count_distribution(selected_counts, selected_metadata, min_counts_per_event, "{}")
	```'''.format(column)))

def all2all():

	return(dedent(
	'''
	```{r, reproducibility, eval=include_all2all, echo=FALSE, results='asis', warning=FALSE}
	cat('### Reproducibility\\n')
	

	if (ncol(selected_counts) <= 10) {
	  cat('To check the reproducibility of biological replicates, we use all2all plots.\\n')
	  all2all(selected_counts_filtered, cex=1)
	} else {
	  cat('With more than 10 samples, meaningful reproducibility plots cannot be created due to size restrictions.')
	}
	```'''))

def run_pca(args):

	return(dedent(
	'''
	```{{r, pca, eval=include_pca, echo=FALSE, results='asis', warning=FALSE}}
	transformation = '{}'
	
	cat('### PCA {{ .tabset .tabset-pills }}\\n')
	if (include_batch_correction) {{
	  cat('#### Pre batch correction\\n')
	}}
	
	pca = run_pca(all_counts_filtered, transformation = transformation)
	print(pca_plot(pca, all_metadata, color_by='{}', shape_by='{}', fill_by='{}', alpha_by='{}', label_by='{}'))
	print(scree_plot(pca))

	if (include_batch_correction) {{
	  cat('\\n\\n#### Post batch correction\\n')
	  pca = run_pca(all_counts_batch_corrected, transformation = transformation)
	  print(pca_plot(pca, all_metadata, color_by='{}', shape_by='{}', fill_by='{}', alpha_by='{}', label_by='{}'))
	  print(scree_plot(pca))
	}}
	```'''.format(args.transformation, args.pca_color, args.pca_shape, args.pca_fill, args.pca_transparency, args.pca_label, args.pca_color, args.pca_shape, args.pca_fill, args.pca_transparency, args.pca_label)))

def run_deseq(treatment, control, highlighted, args, name, design, column):

	return(dedent(
	'''
	```{{r, DESeq_info, eval=include_deseq2, echo=FALSE, results='asis'}}
	cat('## DESeq Analysis {{ .tabset .tabset-pills }}\\n')
	cat('The goal of Differential gene expression analysis is to find genes or transcripts whose difference in expression, when accounting for the variance within condition, is higher than expected by chance.')
	```

	```{{r, run_DESeq, eval=include_deseq2, include=include_deseq2, warning=FALSE, message=FALSE}}
	design = {}
	fitType = "{}"
	padj_significance_cutoff = {}
	fc_significance_cutoff = {}
	padj_floor = {}
	fc_ceiling = {}
	num_labeled = {}
	
	treatment_name = '{}'
	control_name = '{}'
	
	highlighted = c({})

	selected_metadata = selected_metadata %>% mutate({} = factor({}, levels=c(control_name, treatment_name)))

	use_batch_correction_in_DE = {}

	if (include_batch_correction & use_batch_correction_in_DE) {{
	  deseq_input = selected_counts_batch_corrected
	}} else {{
	  deseq_input = selected_counts_filtered 
	}}

	dds = DESeqDataSetFromMatrix(countData = deseq_input, colData = selected_metadata, design = design)
	res <- DESeq(dds, fitType=fitType)
	results = as.data.frame(results(res)) %>% mutate(padj = replace_na(padj, 1)) %>% tibble::rownames_to_column('gene')
	
	if (apply_shrinkage) {{
	  shrinkage_type = "{}"
	  shrink_results = as.data.frame(lfcShrink(res, coef=2, type=shrinkage_type)) %>%
	                   tibble::rownames_to_column('gene') %>%
	                   select(gene, log2FoldChange_shrink=log2FoldChange, lfcSE_shrink=lfcSE)
	  results = results %>% left_join(shrink_results, by='gene')
	}}
	
	post_res = post_processing(results, padj_significance_cutoff=padj_significance_cutoff, fc_significance_cutoff=fc_significance_cutoff, num_labeled=num_labeled, highlighted=highlighted, add_alias={}, fromType='{}',  toType='{}', org='org.Hs.eg.db', apply_shrinkage=apply_shrinkage)
	deseq_results_table(post_res, apply_shrinkage, '{}')
	```'''.format(design, args.fitType, args.padj_significance_cutoff, args.fc_significance_cutoff, args.padj_floor, args.fc_ceiling, args.num_labeled, treatment, control, highlighted, column, column, args.use_bc_in_DE, args.shrinkage_type, args.convert_names, args.count_file_names, args.converted_names, name)))

def ma_plot():
	return(dedent(
	'''
	```{r, ma_plot, echo=FALSE, eval=include_ma & include_deseq2, include=include_ma & include_deseq2, results='asis', warning=FALSE}
	if(include_ma_highlighted) {
	  cat('### MA {.tabset} \\n')
	  cat('#### All Detected \\n')
	} else {
	  cat('### MA\\n')
	}
	
	if(apply_shrinkage) {
	  ma_plot(post_res, Y='log2FoldChange_shrink', padj_cutoff=padj_significance_cutoff, fc_cutoff=fc_significance_cutoff, padj_floor=padj_floor, fc_ceiling=fc_ceiling, num_labeled=num_labeled)
	} else {
	  ma_plot(post_res, padj_cutoff=padj_significance_cutoff, fc_cutoff=fc_significance_cutoff, padj_floor=padj_floor, fc_ceiling=fc_ceiling, num_labeled=num_labeled)
	}

	if(include_ma_highlighted) {
	  cat('\\n\\n#### Highlighted Only \\n')
	  if(apply_shrinkage) {
	    ma_plot(post_res %>% filter(Highlighted==TRUE), Y='log2FoldChange_shrink', padj_cutoff=padj_significance_cutoff, fc_cutoff=fc_significance_cutoff, padj_floor=padj_floor, fc_ceiling=fc_ceiling, num_labeled=num_labeled)
	  } else {
	    ma_plot(post_res %>% filter(Highlighted==TRUE), padj_cutoff=padj_significance_cutoff, fc_cutoff=fc_significance_cutoff, padj_floor=padj_floor, fc_ceiling=fc_ceiling, num_labeled=num_labeled)
	  }
	}
	```'''))

def volcano_plot():
	return(dedent(
	'''
	```{r, volcano_plot, echo=FALSE, eval=include_volcano & include_deseq2, include=include_volcano & include_deseq2, results='asis', warning=FALSE}
	
	if(include_volcano_highlighted) {
	  cat('### Volcano {.tabset} \\n')
	  cat('#### All Detected \\n')
	} else {
	  cat('### Volcano\\n')
	}

	if(apply_shrinkage) {
	  volcano_plot(post_res, X='log2FoldChange_shrink', padj_cutoff=padj_significance_cutoff, fc_cutoff=fc_significance_cutoff, padj_floor=padj_floor, fc_ceiling=fc_ceiling, num_labeled=num_labeled)
	} else {
	  volcano_plot(post_res, padj_cutoff=padj_significance_cutoff, fc_cutoff=fc_significance_cutoff, padj_floor=padj_floor, fc_ceiling=fc_ceiling, num_labeled=num_labeled)
	}

	if(include_volcano_highlighted) {
	  cat('\\n\\n#### Highlighted Only \\n')
	  if(apply_shrinkage) {
	    volcano_plot(post_res %>% filter(Highlighted==TRUE), X='log2FoldChange_shrink', padj_cutoff=padj_significance_cutoff, fc_cutoff=fc_significance_cutoff, padj_floor=padj_floor, fc_ceiling=fc_ceiling, num_labeled=num_labeled)
	  } else {
	    volcano_plot(post_res %>% filter(Highlighted==TRUE), padj_cutoff=padj_significance_cutoff, fc_cutoff=fc_significance_cutoff, padj_floor=padj_floor, fc_ceiling=fc_ceiling, num_labeled=num_labeled)
	  }
	}
	```'''))

def heatmap():
	return(dedent(
	'''
	```{r, heatmap, echo=FALSE, eval=include_heatmap & include_deseq2, include=include_heatmap & include_deseq2, results='asis'}
	cat('### Heatmap\\n')
	
	significant = (post_res %>% filter(Significant == 'Significant'))$gene
	
	heatmap_data = as.data.frame(selected_counts) %>%
	               tibble::rownames_to_column('gene') %>%
	               filter(gene %in% significant) %>%
	               tibble::column_to_rownames('gene')
	
	if (nrow(heatmap_data) > 0) {
	  
	  bad_samples = check_samples(heatmap_data, 1)
	  
	  if (length(bad_samples) == 0) {
	    heatmap_plot(heatmap_data)
	  } else {
	    cat("Cannot include the following samples because they have zero counts across all significant genes which breaks the normalization function:\\n")
	    cat(bad_samples)
	    heatmap_data_filtered = sample_filter(heatmap_data, bad_samples)
	    heatmap_plot(heatmap_data_filtered)
	  }
	} else {
	  cat("No significant genes. Skipping heatmap.")
	}
	```'''))

def write_tables(name):

	return(dedent(
	'''
	```{{r, echo=FALSE, write_tables}}
	prefix = '{}'
	
	if (include_batch_correction) {{
	  write.table(as.data.frame(all_counts) %>% tibble::rownames_to_column('gene'), file=paste0("outputs/", prefix, "_batch_corrected_counts.tsv"), quote=FALSE, row.names=FALSE, sep='\\t')
	}}
	
	if (include_deseq2) {{
	  write.table(results, file=paste0("outputs/", prefix, "_all_deseq2_results.tsv"), quote=FALSE, row.names=FALSE, sep='\\t')
	  write.table(post_res %>% filter(Significant == 'Significant') %>% select(-alias, -Highlighted, -Direction, -Significant, -Label, -Group), file=paste0("outputs/", prefix, "_sig_deseq2_results.tsv"), quote=FALSE, row.names=FALSE, sep='\\t')
	  write.table(post_res, file=paste0("outputs/", prefix, "_results_processed.tsv"), quote=FALSE, row.names=FALSE, sep='\\t')
	}}
	
	if (!include_pca) {{
	  transformation = NA
	}}
	
	if (!include_batch_correction) {{
	  batch_correction_column = NA
	  batch_correction_group_column = NA
	  batch_norm_method = NA
	}}
	
	if (!include_deseq2) {{
	  design = NA
	  fitType = NA
	  padj_significance_cutoff = NA
	  fc_significance_cutoff = NA
	  padj_floor = NA
	  fc_ceiling = NA
	}}
	
	if(!apply_shrinkage) {{
	  shrinkage_type = NA
	}}
	
	output_parameters(prefix, min_counts_per_event, min_samples_per_event, transformation, include_batch_correction, batch_correction_column, batch_correction_group_column, batch_norm_method, design, fitType, apply_shrinkage, shrinkage_type, padj_significance_cutoff, fc_significance_cutoff, padj_floor, fc_ceiling)
	```'''.format(name)))

def write_warnings():

	return(dedent(
	'''
	```{r, write_warnings, echo=FALSE}
	fail_pca = paste("WARNING -", bad_samples, "is not included in PCA because it has less than", min_counts_per_sample, "counts.", sep=' ')
	fail_DE = paste("WARNING -", c(setdiff(treatments, good_treatments), setdiff(controls, good_controls)), "is not included in DE because it has less than", min_counts_per_sample, "counts.", sep=' ')
	writeLines(c(fail_pca, fail_DE), paste0("outputs/", prefix, "_warning_log.txt"))
	```'''))

def session_info():

	return(dedent(
	'''
	```{r, session_info, echo=FALSE, results='asis'}
	cat('## Session Info\\n')
	sessionInfo()
	```'''))

def write_output(name, output):

	with open('%s.Rmd' % (name), 'w') as out:
		out.write('\n'.join(output)) 

def parseArguments():
	parser = argparse.ArgumentParser(prog="prepare_DESeq2", description='', usage='%(prog)s [options]')
	
	input_args = parser.add_argument_group('Input')
	input_args.add_argument('-c', '--counts', required=True, help='Name of count file.', metavar='', dest='counts')
	input_args.add_argument('-g', '--groups', required=True, help='Name of groups file.', metavar='', dest='groups')
	input_args.add_argument('-x', '--comparisons', required=True, help='Name of comparisons file.', metavar='', dest='comparisons')
	input_args.add_argument('-p', '--prefix', default='', help='Path to input directory.', metavar='', dest='input_prefix')


	qc_modules = parser.add_argument_group('QC Modules')
	qc_modules.add_argument('--include-distribution', choices=['TRUE', "FALSE"], default='TRUE', help='Create count distribution plot', metavar='', dest='include_distribution')
	qc_modules.add_argument('--include-all2all', choices=['TRUE', "FALSE"], default='TRUE', help='Create reproducibility plots', metavar='', dest='include_all2all')
	qc_modules.add_argument('--include-pca', choices=['TRUE', "FALSE"], default='TRUE', help='Create PCA plots', metavar='', dest='include_pca')

	filtering = parser.add_argument_group('Filtering Criteria')
	filtering.add_argument('--filter-type', choices=["global","local"], default='local', help='Apply filter across all samples (global) or just samples within the specific comparison (local).', metavar='', dest='filter_type')
	filtering.add_argument('--min-counts-per-event', type=int, default=10, help='Filter genes with less than this count.', metavar='', dest='min_counts_per_event')
	filtering.add_argument('--min-samples-per-event', type=int, default=1, help='Filter genes which do not have at least min-counts-per-event reads in this many samples.', metavar='', dest='min_samples_per_event')
	filtering.add_argument('--min-counts-per-sample', type=int, default=1000, help='Filter samples which do not have at least min_counts_per_sample reads.', metavar='', dest='min_counts_per_sample')
	filtering.add_argument('--excluded-events', type=str, nargs='+', default='', help='List of events to remove from analysis.', metavar='', dest='excluded_events')

	pca = parser.add_argument_group('PCA Settings')
	pca.add_argument('--transformation', choices=["Default","rlog","vst","None"], default='Default', help='Transformation to apply prior to PCA (does not apply to data sent into DESeq). Default will use the rlog transformation for up to 50 samples and vst if there are more than 50 samples.', metavar='', dest='transformation')
	pca.add_argument('--pca-color', type=str, default='', help='Column from metadata file by which to determine point color.', metavar='', dest='pca_color')
	pca.add_argument('--pca-shape', type=str, default='', help='Column from metadata file by which to determine point shape.', metavar='', dest='pca_shape')
	pca.add_argument('--pca-fill', type=str, default='', help='Column from metadata file by which to determine point fill.', metavar='', dest='pca_fill')
	pca.add_argument('--pca-transparency', type=str, default='', help='Column from metadata file by which to determine point transparency.', metavar='', dest='pca_transparency')
	pca.add_argument('--pca-label', type=str, default='', help='Column from metadata file by which to determine point label.', metavar='', dest='pca_label')

	batch_correction = parser.add_argument_group('Batch Correction Settings')
	batch_correction.add_argument('--include-batch-correction', choices=["TRUE", "FALSE"], default='FALSE', help='Perform batch correction prior to DESeq', metavar='', dest='include_batch_correction')
	batch_correction.add_argument('--batch-correction-column', type=str, default='batch', help='Column used for batch correction.', metavar='', dest='batch_correction_column')
	batch_correction.add_argument('--batch-correction-group-column', type=str, default='group', help='Column containing biological variable grouping.', metavar='', dest='batch_correction_group_column')
	batch_correction.add_argument('--batch-normalization-algorithm', choices=["MRN","TMM","RLE","upperquartile", "none"], default='TMM', help='Normalization algorithm to use prior to batch correction.', metavar='', dest='batch_correction_normalization_algorithm')

	deseq2 = parser.add_argument_group('DESeq2 Settings')
	deseq2.add_argument('--include-DESeq2', choices=["TRUE", "FALSE"], default='TRUE', help='Run DESeq2', metavar='', dest='include_DESeq2')
	deseq2.add_argument('--design', type=str, default='~group', help='Design formula', metavar='', dest='design')
	deseq2.add_argument('--fitType', choices=["parametric","local","mean","glmGamPoi"], default='parametric', help='fitType used during DESeq2 analysis.', metavar='', dest='fitType')
	deseq2.add_argument('--apply-shrinkage', choices=["TRUE", "FALSE"], default='TRUE', help='Apply shrinkage to fold change after DESeq2', metavar='', dest='apply_shrinkage')
	deseq2.add_argument('--use-batch-correction-in-DE', choices=["TRUE", "FALSE"], default='TRUE', help='Use batch corrected values in DESeq', metavar='', dest='use_bc_in_DE')
	deseq2.add_argument('--shrinkage-type', choices=["apeglm","ashr","normal"], default='apeglm', help='Shinkage algorithm to use.', metavar='', dest='shrinkage_type')
	deseq2.add_argument('--include-volcano', choices=["TRUE", "FALSE"], default='TRUE', help='Create volcano plots', metavar='', dest='include_volcano')
	deseq2.add_argument('--include-ma', choices=["TRUE", "FALSE"], default='TRUE', help='Create MA plots', metavar='', dest='include_ma')
	deseq2.add_argument('--include-heatmap', choices=["TRUE", "FALSE"], default='TRUE', help='Create heatmap', metavar='', dest='include_heatmap')

	significance = parser.add_argument_group('Significance Settings')
	significance.add_argument('--padj-significance-cutoff', type=float, default=.05, help='Adjusted p-value significance threshold required to be called differentially expressed.', metavar='', dest='padj_significance_cutoff')
	significance.add_argument('--fc-significance-cutoff', type=float, default=1, help='log2 Fold Change threshold required to be called differentially expressed.', metavar='', dest='fc_significance_cutoff')
	significance.add_argument('--padj-floor', type=str, default='0', help='Reassign adjusted p-values below this value to this value. Set at 0 to disable.', metavar='', dest='padj_floor')
	significance.add_argument('--fc-ceiling', type=str, default='Inf', help='Reassign fold changes (log2) greater than this value to this value. Set to Inf to disable.', metavar='', dest='fc_ceiling')

	label = parser.add_argument_group('Label Options')
	label.add_argument('--convert-names', choices=['TRUE', 'FALSE'], default='FALSE', help='Convert gene names', metavar='', dest='convert_names')
	label.add_argument('--count-file-names', choices=["ENSEMBL","ENTREZID","SYMBOL"], default='SYMBOL', help='Name used in count file.', metavar='', dest='count_file_names')
	label.add_argument('--converted-names', choices=["ENSEMBL","ENTREZID","SYMBOL"], default='ENSEMBL', help='Name to use in plots.', metavar='', dest='converted_names')
	label.add_argument('--num-labeled', type=str, default='Inf', help='Number of significant genes to label. Set to Inf to label all. Set to 0 to turn labels off.', metavar='', dest='num_labeled')
	label.add_argument('--highlighted-genes', type=str, nargs='+', default='', help='Space separated list of genes to label regardless of significance.', metavar='', dest='highlighted_genes')
	label.add_argument('--include-volcano-highlighted', choices=["TRUE", "FALSE"], default='FALSE', help='Create volcano plots of highlighted genes only.', metavar='', dest='include_volcano_highlighted')
	label.add_argument('--include-ma-highlighted', choices=["TRUE", "FALSE"], default='FALSE', help='Create MA plots of highlighted genes only.', metavar='', dest='include_ma_highlighted')

	return parser.parse_args()

if __name__ == "__main__":
	args = parseArguments()
	main(args)
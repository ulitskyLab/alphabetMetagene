#!/usr/bin/env Rscript


library(ggplot2)
library(optparse)
library(reshape2)
library(stringr)
library(RColorBrewer)


EXTENDED_BED_FILE_SIZE = 7
NAME_MATCH_COLUMNS_SIZE = 2
ACCURACY_FREQ = 3
ACCURACY_SD = 2
ACCURACY_PVAL = 8
SIGNIFICANCE_LEVELS = c(0.0001, 0.001, 0.01, 0.05, 0.1, 1)


option_list = list(
  make_option(c("-f", "--fasta_paths"), type="character", default=NULL, 
              help="fasta paths vector (sep by commas), \
              first file should be the sample file to which others are compared", metavar="character"),
  make_option(c("-m", "--matching_control_files"), type="character", default=NULL,
              help="bed file per control with a 7th column containing the name of the matching sample entry.\
              Alternatively, a file that contains the name in the first column and what it controls for in the second (separated by tabs)", metavar="character"),
  make_option(c("-a", "--alphabet_vec"), type="character", default='DNA',
              help="alphabet used (either DNA/RNA/protein or alphabet string separated by commas) \
              [default= %default]", metavar="character"),
  make_option(c("-t", "--title"), type="character", default='',
              help="output plot title (defaults is none)", metavar="character"),
  make_option(c("-s", "--sample_names"), type="character", default=NULL,
              help="sample and control names (defaults to their index numbers)", metavar="character"),
  make_option(c("-w", "--plot_width"), type="numeric", default=12,
              help="plot width [default= %default] (relevant for pdf default output option)", metavar="numeric"),
  make_option(c("-l", "--plot_height"), type="numeric", default=NULL,
              help="plot height [default= %default]", metavar="numeric"),
  make_option(c("-c", "--colors"), type="character", default=NULL,
              help="plot colors, \
              can use color names (e.g. 'black') \
              or color coding (e.g. = #FF0000) \
              separated by commas \
              at the length of the alphabet", metavar="character"),
  make_option(c("-x", "--text_size"), type="numeric", default=16,
              help="bars text size", metavar="numeric"),
  make_option(c("-u", "--make_upper"), type="logical", default=TRUE,
              help="whether the fasta files should be moved to capital letters", metavar="logical"),
  make_option(c("-r", "--format_output"), type="character", default='pdf',
              help="format of the output file (available are pdf/png/jpeg)", metavar="character"),
  make_option(c("-o", "--out_path"), type="character", default='out.pdf', 
              help="output file name [default= %default] defaults to current directory", metavar="character")
  ); 



opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);



# assumption: first file is the sample file, following are the control files
read_fastas = function(fasta_path_vec, make_upper) {
  fasta_file_path = unlist(strsplit(fasta_path_vec, split = ","))

  fasta_files = list()
  
  for (i in 1:length(fasta_file_path)) {
    
    fasta_file = read.csv(fasta_file_path[i], stringsAsFactors = FALSE, header = FALSE)
    fasta_files[[i]] = as.data.frame(gsub(" .*", "", gsub(">", "", fasta_file[c(TRUE, FALSE),])))
    fasta_files[[i]]$sequence = fasta_file[c(FALSE, TRUE),]
    if (make_upper) {
      fasta_files[[i]]$sequence = toupper(fasta_files[[i]]$sequence)
    }
    colnames(fasta_files[[i]]) = c('name', 'sequence')
    
  }
  return(fasta_files)
}



# in cases of matching controls per sample
read_matching_control_data = function(matching_sample_control_path_vec) {
  matching_samples_path = unlist(strsplit(matching_sample_control_path_vec, split = ","))
  
  matching_files = list()
  
  for (i in 1:length(matching_samples_path)) {
        
    matching_files[[i]] = read.table(matching_samples_path[i], stringsAsFactors = FALSE, header = FALSE)
    
    if (ncol(matching_files[[i]]) == EXTENDED_BED_FILE_SIZE) {
      colnames(matching_files[[i]]) = c('chr', 'start', 'end', 'name', 'score', 'strand', 'control_for')
    }
    
    else if (ncol(matching_files[[i]]) == NAME_MATCH_COLUMNS_SIZE) {
      colnames(matching_files[[i]]) = c('name', 'control_for')
    }
    
    else {
      print('Input match control files are not in the right dimensions (should be either 2 or 7 columns).')
      return(FALSE)
    }
  }
  return(matching_files)
}



# alphabet vec is a character vector of single letters of the alphabet
calc_alphabet_ratio = function(fastas_df_list, alphabet_vec) {

  for (i in 1:length(fastas_df_list)) {
    
    fastas_df_list[[i]]$seq_length = nchar(fastas_df_list[[i]]$sequence)
    fastas_df_list[[i]][alphabet_vec] = 0

    for (j in 1:length(alphabet_vec)) {
      fastas_df_list[[i]][,alphabet_vec[j]] = str_count(fastas_df_list[[i]]$sequence, alphabet_vec[j]) / fastas_df_list[[i]][,'seq_length']
      # fastas_df_list[[i]][,alphabet_vec[j]] = lengths(regmatches(fastas_df_list[[i]]$sequence, gregexpr(as.character(alphabet_vec[j]), fastas_df_list[[i]]$sequence))) / fastas_df_list[[i]][,'seq_length']
    }
  }
  return(fastas_df_list)
}



aggregate_control_ratios = function(control_dfs, matching_sample_files, alphabet_vec) {
  
  aggregated_controls = list()
  
  for (i in 1:length(control_dfs)) {
 
   # match
    control_dfs[[i]]$control_for = apply(control_dfs[[i]],
                                         1,
                                         function(x) {matching_sample_files[[i]][matching_sample_files[[i]]$name == x['name'], 'control_for']})    
    
    # split
    ctrl_sub_mat_split = split(control_dfs[[i]], control_dfs[[i]]$control_for)

    # aggregate
    chunks_agg = lapply(ctrl_sub_mat_split, function(x) {sapply(x[,alphabet_vec], mean)})

    # convert back to df
    chunks_agg_united = as.data.frame(do.call(rbind, chunks_agg))
    chunks_agg_united$matched_sample = row.names(chunks_agg_united)
    row.names(chunks_agg_united) = NULL
    
    aggregated_controls[[i]] = chunks_agg_united

  }
  return(aggregated_controls)
}



calc_wilcoxon = function(agg_controls_dfs, sample_df, alphabet_vec, matching=T) {
  sample_df = sample_df[,c('name', alphabet_vec)]
  pvals = sample_df[F, alphabet_vec]

  if (matching) {

    for (i in 1:length(agg_controls_dfs)) {
      agg_controls_dfs[[i]] = agg_controls_dfs[[i]][order(match(agg_controls_dfs[[i]]$matched_sample, sample_df$name)),]
      for (j in 1:length(alphabet_vec)) {
        pvals[i, alphabet_vec[j]] = wilcox.test(sample_df[,alphabet_vec[j]], agg_controls_dfs[[i]][,alphabet_vec[j]], paired=TRUE)$p.value
      }
    }
  }
  
  else { # in that case the control_dfs will actually not be aggregated..
    for (i in 1:length(agg_controls_dfs)) {
      agg_controls_dfs[[i]] = agg_controls_dfs[[i]][,alphabet_vec]
      for (j in 1:length(alphabet_vec)) {
        pvals[i, alphabet_vec[j]] = wilcox.test(sample_df[,alphabet_vec[j]], agg_controls_dfs[[i]][,alphabet_vec[j]])$p.value
      }
    }
  }
  pvals$compared_to_origin = as.numeric(rownames(pvals)) + 1
  return(pvals)
}



aggregate_frequencies = function(fastas_dfs, alphabet_vec) {
  frequencies = fastas_dfs[[1]][F, alphabet_vec]
  sds = frequencies
  
  for (i in 1:length(fastas_dfs)) {
    frequencies[i,alphabet_vec] = apply(fastas_dfs[[i]][,alphabet_vec], 2, mean)
    sds[i,alphabet_vec] = apply(fastas_dfs[[i]][,alphabet_vec], 2, sd)
    frequencies[i, 'origin'] = as.character(i)
    sds[i, 'origin'] = as.character(i)
  }
  
  return(list(frequencies, sds))
}



plot_frequencies = function(frequencies_df, sds, pvals, title, colors, sample_names=unique(frequencies_df$origin), text_size) {
  
  frequencies_df = melt(frequencies_df)
  num_cat = length(unique(frequencies_df$variable))
  sds = melt(sds)
  frequencies_df$round_val = round(frequencies_df$value, ACCURACY_FREQ)
  freq = merge(frequencies_df, sds, by=c('origin', 'variable'))
  freq$round_sds = round(freq$value.y, ACCURACY_SD)

  if (!is.null(text_size)) {
    text_size_final = text_size
  }
  else {
    text_size_final = (text_size/(num_cat**2)) + 2
  }
    
  p = ggplot(freq, aes(x=origin, y=value.x, fill=variable, label=round_val)) + 
    geom_bar(stat="identity") + 
    geom_text(size = text_size_final+1, position = position_stack(vjust = 0.5)) + 
    geom_text(aes(label=paste0('SD=',round_sds)), size=text_size_final, position = position_stack(vjust = 0.35)) +
    ylab("Frequency") + #scale_x_discrete(labels=c("1" = "tandem", "2" = "promoter", "3" = "3prime")) +
    ggtitle(title) + scale_x_discrete(labels=sample_names) + theme_bw() + scale_fill_manual(values=colors)
  
  if (length(sample_names) > 1) {
    pvals_df2 = melt(pvals, id.vars = c('compared_to_origin'))
    
    freq = merge(freq, pvals_df2, by.x=c('origin', 'variable'), by.y=c('compared_to_origin', 'variable'), all=TRUE)
    freq$round_pval = round(freq$value, ACCURACY_PVAL)
    freq$sig_level = ifelse(freq$value > SIGNIFICANCE_LEVELS[4], 'ns',
                                ifelse(((freq$value <= SIGNIFICANCE_LEVELS[4]) & (freq$value > SIGNIFICANCE_LEVELS[3])), '*', 
                                       ifelse(((freq$value <= SIGNIFICANCE_LEVELS[3]) & (freq$value > SIGNIFICANCE_LEVELS[2])), '**', 
                                              ifelse(((freq$value <= SIGNIFICANCE_LEVELS[2]) & (freq$value > SIGNIFICANCE_LEVELS[1])), '***', '****'))))
    
    p = p + geom_text(data=freq, aes(label=paste0('P=',sig_level)), size=text_size_final, position = position_stack(vjust = 0.2))
  }
  return(p)
}



save_plot = function(p, out_path, width, height, format_output) {
  if (format_output == 'jpeg') {
    jpeg(out_path)
  }
  else if (format_output == 'pdf') {
    pdf(out_path, width)
  }
  else {
    png(out_path)
  }
  print(p)
  dev.off()
  return()
}



parse_input_alphabet = function(alphabet_vec) {
  
  if (alphabet_vec == 'DNA') {
    alphabet_vec = c('A','C','G','T')
  }
  else if (alphabet_vec == 'RNA') {
    alphabet_vec = c('A','C','G','U')
  }
  else if (alphabet_vec == 'protein') {
    alphabet_vec = c('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y')
  }
  else {
    alphabet_vec = unlist(strsplit(alphabet_vec, split = ","))
  }
  return(alphabet_vec)
}



overall_frequency = function(fastas_vec,
                             matching_controls_vec,
                             alphabet_vec,
                             title,
                             sample_names,
                             plot_width,
                             plot_height,
                             colors,
                             text_size,
                             make_upper,
                             format_output,
                             out_path) {
  
  alphabet_vec = parse_input_alphabet(alphabet_vec)
  fastas = read_fastas(fastas_vec, make_upper)
  
  alphabet_ratios = calc_alphabet_ratio(fastas, alphabet_vec)

  controls_fastas = alphabet_ratios[2:length(alphabet_ratios)]
  sample_fastas = alphabet_ratios[[1]]
  
  if (!(is.null(matching_controls_vec))) {

    matching_controls = read_matching_control_data(matching_controls_vec)
    agg_controls = aggregate_control_ratios(controls_fastas, matching_controls, alphabet_vec)
    sample_df = list(sample_fastas[,c(alphabet_vec, 'name')])
    ratios_agg_df = c(sample_df, agg_controls)    
    frequencies_data = aggregate_frequencies(ratios_agg_df, alphabet_vec)
    pvals_df = calc_wilcoxon(agg_controls, sample_fastas, alphabet_vec)
  }

  
  else if ((length(fastas) > 1)) {
    pvals_df = calc_wilcoxon(controls_fastas, sample_fastas, alphabet_vec, F)
    frequencies_data = aggregate_frequencies(alphabet_ratios, alphabet_vec)
  }

  
  # in case there's only one fasta file as input no test is done
  else {
    pvals_df = NULL
    frequencies_data = aggregate_frequencies(alphabet_ratios, alphabet_vec)
  }
  
  frequencies = frequencies_data[[1]]
  sds = frequencies_data[[2]]
  
  if (!(is.null(sample_names))) {
    sample_names = unlist(strsplit(sample_names, split = ","))
  }
  else {
    sample_names = paste0(rep('sample', length(fastas)), as.character(seq(1,length(fastas),1)))
  }  
  
  if (!is.null(colors)) {
    color_names = unlist(strsplit(colors, split = ","))
  }
  else {
    coul = brewer.pal(4, "PuOr") 
    color_names = colorRampPalette(coul)(length(alphabet_vec))
  }
  
  p = plot_frequencies(frequencies, sds, pvals_df, title, color_names, sample_names, text_size)
  save_plot(p, out_path, plot_width, plot_height, format_output)
  
  return(print("finished succesfully"))
}



overall_frequency(opt$fasta_paths,
                  opt$matching_control_files,
                  opt$alphabet_vec,
                  opt$title,
                  opt$sample_names,
                  opt$plot_width,
                  opt$plot_height,
                  opt$colors,
                  opt$text_size,
                  opt$make_upper,
                  opt$format_output,
                  opt$out_path)

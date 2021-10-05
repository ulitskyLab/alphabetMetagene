#!/usr/bin/env Rscript

library(optparse)
library(plyr)
library(stringr)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(pheatmap)
library(gridExtra)


ALTERNATIVE_NUMBER_BINS = 100
SIGNIFICANCE_LEVELS = c(0, 0.001, 0.01, 0.05, 0.1, 1)
EXTENDED_BED_FILE_SIZE = 7
NAME_MATCH_COLUMNS_SIZE = 2

FIRST_BIN_SIZE_COL_INDEX = 8
MIDDLE_BIN_SIZE_COL_INDEX = 5
LAST_BIN_SIZE_COL_INDEX = 9
START_OF_BINS_SEQUENCE_INDEX = 10

ZERO_TO_ONE_PVAL_SCALE = c('black','gray15','gray58', 'gray79','gray90','gray95','gray98','gray99','white')
PVAL_COLORS = c('black','gray20','gray50','gray70','white')


option_list = list(
  make_option(c("-f", "--fasta_path_vec"), type="character", default=NULL, 
              help="fasta paths vec (sep by commas)", metavar="character"),
  make_option(c("-m", "--match"), type="character", default=NULL,
              help="Input options: \
              1. Ignore -m option [default= %default]: will calculate non-paired Wilcoxon rank-sum test between samples. \
              2. Input .bed file per control with an extra 7th column containing the name of the matching sample entry.\
              3. Alternatively, a file that contains the name in the first column and what it controls for in the second (separated by tabs)", metavar="character"),
  make_option(c("-a", "--alphabet_vec"), type="character", default='A,C,G,T', 
              help="alphabet used (either DNA/RNA/protein or alphabet string separated by commas) [default= %default]", metavar="character"),
  make_option(c("-s", "--sample_names_vec"), type="character", default=NULL,
              help="sample and control names (defaults to their index numbers)", metavar="character"),
  make_option(c("-t", "--output_title"), type="character", default=NULL,
              help="output plot title (Prefix is the alphabet character)", metavar="character"),
  make_option(c("-o", "--out_full_path"), type="character", default=NULL, 
              help="output file name (defaults to current directory)", metavar="character"),
  make_option(c("-u", "--make_upper"), type="logical", default=TRUE,
              help="whether the fasta files should be moved to capital letters", metavar="logical"),
  make_option(c("-c", "--plot_colors"), type="character", default=NULL, 
              help="color names separated by comma (should be as the length of sample names)", metavar="character")
); 



opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);



# assumption: first file is the sample file, following are the control files
read_fastas = function(fasta_path_vec, make_upper=TRUE) {
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
      print('Input match control files are not in the right dimensions (should have either 2 or 7 columns).')
      return(FALSE)
    }
  }
  return(matching_files)
}



bind_fastas = function(fasta_files_lst) {
  fastas = mapply(cbind, fasta_files_lst, "origin"=as.character(1:length(fasta_files_lst)), SIMPLIFY = F)
  bound_fastas = ldply(fastas, rbind)
  return(bound_fastas)
}



define_num_bins = function(combined_fasta_file) {
  
  combined_fasta_file$seq_length = as.numeric(sapply(combined_fasta_file$sequence, nchar))
  
  if (length(unique(combined_fasta_file$seq_length)) == 1) {
    num_bins = combined_fasta_file$seq_length[1]
  }
  
  else {
    num_bins = min(combined_fasta_file$seq_length, ALTERNATIVE_NUMBER_BINS)
  }
  return(list(num_bins, combined_fasta_file))
}



get_num_bases_per_bin = function(combined_fasta, num_bins) {  
  combined_fasta$middle_bin_size = floor(combined_fasta$seq_length / num_bins)
  combined_fasta$total_middle_bins_length = combined_fasta$middle_bin_size * (num_bins - 2)
  combined_fasta$first_last_bin_sizes = ((combined_fasta$seq_length - combined_fasta$total_middle_bins_length) / 2)
  return(combined_fasta)
}



divide_bins = function(combined_fasta) {
  # calc exact first and last bin sizes
  combined_fasta$first_bin_sizes = floor(combined_fasta$first_last_bin_sizes)
  combined_fasta$last_bin_sizes = ceiling(combined_fasta$first_last_bin_sizes)
  
  combined_fasta$first_bin = mapply(function(x,y) {substr(x,1,y)}, x=combined_fasta$sequence, y=combined_fasta$first_bin_sizes)
  combined_fasta$last_bin = mapply(function(x,y) {substr(x,nchar(x)-y+1,nchar(x))}, x=combined_fasta$sequence, y=combined_fasta$first_bin_sizes)
  combined_fasta$all_middle_bins = mapply(function(x,y,z) {substr(x,y+1,nchar(x)-z)}, x=combined_fasta$sequence, y=combined_fasta$first_bin_sizes, z=combined_fasta$last_bin_sizes)
  
  # returns mat:
  start_pos_middle_bins = mapply(function(x,y) {seq(1,nchar(x),by=y)}, x=combined_fasta$all_middle_bins, y=combined_fasta$middle_bin_size)
  # convert to list:
  start_pos_middle_bin_list = split(start_pos_middle_bins, rep(1:ncol(start_pos_middle_bins), each = nrow(start_pos_middle_bins)))
  
  # rows are bin numbers, columns are partitions of specific FASTA entry:
  middle_bins_mat = mapply(function(X,Y,Z) {sapply(Y, function(ii) {substr(X, ii, ii+Z-1)})},
                           X=combined_fasta$all_middle_bins, Y=start_pos_middle_bin_list, Z=combined_fasta$middle_bin_size)
  middle_bins_mat_t = t(middle_bins_mat)
  middle_bins_mat_t = as.data.frame(middle_bins_mat_t)
  colnames(middle_bins_mat_t) = sapply(colnames(middle_bins_mat_t), function(x) {paste0('bin_', x)})
  
  combined_fasta = cbind(combined_fasta, middle_bins_mat_t)
  rownames(combined_fasta) = NULL
  return(combined_fasta)
}



# calculate wilcoxon between sample and controls
calculate_wilcoxon = function(alphabet, combined, sample_names, matching=NULL, num_bin=NULL) {

  rownames(combined) = paste(combined$name, combined$origin, sep="_")

  combined[,c(1,3,13:(length(colnames(combined))))] = apply(combined[,c(1,3,13:(length(colnames(combined))))], 2, function(x) {as.character(x)})
  
  combined = combined[,c(1:10, 13:length(colnames(combined)), 11)]
  
  # keep: name | origin | bin1 |... | binN
  combined = combined[,c(1,3,10:length(colnames(combined)))]
  
  pvals_df_list = list()

  if (is.null(matching)) {
    for (i in 1:length(alphabet)) {
      # start a pvals df for the current alphabet letter
      
      pvals_df = data.frame(matrix(NA, nrow = (length(sample_names)-1), ncol = (ncol(combined)-2)))

      rownames(pvals_df) = sample_names[2:length(sample_names)]
      colnames(pvals_df) = as.character(seq(1,ncol(pvals_df),1))

      cur_alphabet = combined
      cur_alphabet[,3:ncol(cur_alphabet)][] = lapply(cur_alphabet[,3:ncol(cur_alphabet)], function(x) {str_count(x, alphabet[i])})
      cur_alphabet_per_sample = split(cur_alphabet, cur_alphabet$origin)
      
      for (j in 2:length(sample_names)) {
        for (u in 3:ncol(cur_alphabet)) {
          result = lapply(cur_alphabet_per_sample, "[", , u)
          pvals_df[sample_names[j], u-2] = wilcox.test(result[[1]], result[[j]])$p.value
        }
      }
      pvals_df_list[[i]] = pvals_df
    }
  }

  else {
    # create a df to contain the pvalues
    pvals_df = data.frame(matrix(NA, nrow = (length(sample_names)-1), ncol = (ncol(combined)-2)))

    rownames(pvals_df) = sample_names[2:length(sample_names)]
    colnames(pvals_df) = as.character(seq(1,ncol(pvals_df),1))
    
    colnames(combined)[3:ncol(combined)] = as.character(seq(1,(ncol(combined)-2),1))
    combined$control_for = apply(combined, 1, function(x) {ifelse(as.numeric(x['origin']) == 1,
                                                                              x['name'],
                                                                              matching[[as.numeric(x['origin'])-1]][matching[[as.numeric(x['origin'])-1]]$name == x['name'], 'control_for'])})
    
    # divide based on both origin and what it controls for
    split_var = list(combined$origin, combined$control_for)
    combined_split = split(combined, split_var)


    for (i in 1:length(alphabet)) {
      combined_split_1 = lapply(combined_split, function(x) {apply(x[3:(ncol(x)-1)], 2, function(y) {mean(str_count(y, alphabet[i]))})})
      combine_by_origin = list()

      for (j in 1:length(sample_names)) {
        combine_by_origin[[j]] = combined_split_1[startsWith(names(combined_split_1), as.character(j))]
      }

      dfs = lapply(combine_by_origin, function(x) {t(as.data.frame(x))})
      
      for (j in 2:length(sample_names)) {
        for (u in 1:ncol(dfs[[1]])) {
          result = lapply(dfs, "[", , u)
          pvals_df[sample_names[j], u] = wilcox.test(result[[1]], result[[j]], paired = T)$p.value
        }
      }

      pvals_df_list[[i]] = pvals_df
    }
  }

  if (!is.null(num_bin)) {
    for (i in 1:length(pvals_df_list)) {
      pvals_df_list[[i]][] = lapply(pvals_df_list[[i]], function(x) {p.adjust(as.numeric(x), method = "bonferroni", n=num_bin)})
    }
  }
  return(pvals_df_list)
}



calculate_bin_frequency = function(alphabet, combined, matching=NULL) {

  rownames(combined) = paste(combined$name, combined$origin, sep="_")
  combined[,c(1,3,13:(length(colnames(combined))))] = apply(combined[,c(1,3,13:(length(colnames(combined))))], 2, function(x) {as.character(x)})
  combined = combined[,c(1:10, 13:length(colnames(combined)), 11)]
  ctrl_sub_mat_split_by_origin = split(combined, combined$origin)

  combined_bins = lapply(ctrl_sub_mat_split_by_origin, function(x) {x[,10:(length(colnames(x)))]})

  fasta_nuc = as.data.frame(combined_bins[[1]][F,], stringsAsFactors = F)
  fastas_nuc = rep(list(fasta_nuc), length(combined_bins))

  
  if (is.null(matching)) {
    for (i in 1:length(combined_bins)) {
      for (j in 1:length(alphabet)) {
        fastas_nuc[[i]][alphabet[j],] = apply(combined_bins[[i]],
                                              2,
                                              function(x) {(as.numeric(str_count(paste(x, collapse = ''),
                                                                                        alphabet[j])))})        
      }

      fastas_nuc[[i]] = apply(fastas_nuc[[i]], 2, function(x) {as.numeric(x)/sum(as.numeric(x))})
      rownames(fastas_nuc[[i]]) = alphabet
    }

    fastas_nuc = lapply(fastas_nuc, as.data.frame)
  }

  
  else {
    comb_split = list()

    for (i in 1:length(ctrl_sub_mat_split_by_origin)) {
      ctrl_sub_mat_split_by_origin[[i]]$control_for = apply(ctrl_sub_mat_split_by_origin[[i]], 1, 
                                   function(x) {ifelse(as.numeric(x['origin']) == 1,
                                                       x['name'],
                                                       matching[[as.numeric(x['origin'])-1]][matching[[as.numeric(x['origin'])-1]]$name == x['name'], 'control_for'])})
      
      split_combined_by_control = split(ctrl_sub_mat_split_by_origin[[i]], ctrl_sub_mat_split_by_origin[[i]]$control_for)
      comb_split[[i]] = split_combined_by_control
      
      # aggregate
      for (j in 1:length(alphabet)) {
        chunks_agg = lapply(comb_split[[i]], function(x) {apply(x[,START_OF_BINS_SEQUENCE_INDEX:(length(colnames(combined)))], 2, function(y) {(as.numeric(str_count(paste(y, collapse = ''), alphabet[j])))})})
        # size of first and last bins and size of middle bins (for each group of controls):
        sum_of_chunks = lapply(comb_split[[i]], function(x) {apply(x[,c(FIRST_BIN_SIZE_COL_INDEX,rep(MIDDLE_BIN_SIZE_COL_INDEX,(length(START_OF_BINS_SEQUENCE_INDEX:(length(colnames(combined)))) - 2)),LAST_BIN_SIZE_COL_INDEX)], 2, function(y) {sum(y)})})
        combined_aggregation = as.data.frame(chunks_agg) / as.data.frame(sum_of_chunks)
        chunks_agg2 = apply(combined_aggregation, 1, function(x) {mean(x)})
        fastas_nuc[[i]][alphabet[j],] = chunks_agg2
      }
    }
  }
  return(fastas_nuc)
}



arrange_df_to_plot = function(bin_freq_lst, alphabet, sample_names, pval_dfs) {
  
  all_res = list()

  for (i in 1:length(alphabet)) {
    result = lapply(bin_freq_lst, "[", alphabet[i], )
    results_df = as.data.frame(do.call(rbind, result))

    row.names(results_df) = sample_names
    colnames(results_df) = seq(1, ncol(results_df), 1)
        
    if (length(sample_names) > 1) {
      results_df[paste0(sample_names[1], '.p'),] = 1
      for (j in 2:length(sample_names)) {
        results_df[paste0(sample_names[j], '.p'),] = pval_dfs[[i]][sample_names[j],]
      }
    }
    all_res[[i]] = results_df
  }
  return(all_res)
}



melt_df = function(df_lst, sample_names) {
  
  transformed_df_lst = lapply(df_lst, function(x) {t(x)})
  
  freq_t_df_lst = lapply(transformed_df_lst, "[", , sample_names)
  melt_t = (lapply(freq_t_df_lst, function(x) {melt(x)}))

  melt_merge = list()
  
  if (length(sample_names) > 1) {
    pval_t_df_lst = lapply(transformed_df_lst, "[", , seq(length(sample_names)+1,length(sample_names)*2,1))
    melt_t_pval = (lapply(pval_t_df_lst, function(x) {melt(x)}))

    for (i in 1:length(melt_t_pval)) {

      colnames(melt_t_pval[[i]]) = c('bin', 'p', 'pval_start')
      colnames(melt_t[[i]]) = c('bin', 'sample', 'frequency')

      melt_t_pval[[i]]$sample = apply(melt_t_pval[[i]], 1, function(x) {unlist(strsplit(x[2], '\\.'))[1]})
      melt_merge[[i]] = merge(x=melt_t[[i]], y=melt_t_pval[[i]], by=c('sample', 'bin'))
      melt_merge[[i]] = melt_merge[[i]][order(melt_merge[[i]]$bin), ]
    }
  }

  else {
    for (i in 1:length(melt_t)) {
      colnames(melt_t[[i]]) = c('frequency')
      melt_t[[i]]$bin = rownames(melt_t[[i]])
      melt_t[[i]]$sample = sample_names[1]
      melt_t[[i]] = melt_t[[i]][,c('bin', 'sample', 'frequency')]
      melt_merge[[i]] = melt_t[[i]]
    }
  }
  return(melt_merge)
}


plot_metagene = function(df, alphabet, sample_names, added_title, out_full_path, plot_colors) {

  plot_colors = parse_plot_colors(plot_colors, sample_names)
  
  pdf(out_full_path, width=12)
  

  for (i in 1:length(alphabet)) {

    cur_df = df[[i]]
    cur_df[,'bin'] = as.numeric(as.character(cur_df[,'bin']))
    cur_df[,'frequency'] = as.numeric(as.character(cur_df[,'frequency']))
    
    p = ggplot(cur_df, mapping=aes(x=bin, y=frequency)) +
      geom_line(aes(colour=sample)) +
      geom_point(aes(x=bin, y=frequency, colour=sample), size=1, alpha=0.9) + 
      scale_color_manual(values=plot_colors) +
      labs(y='alphabet ratio', x='bin') +
      ggtitle(paste0(alphabet[i], ' ', added_title)) +
      theme_bw()
    
    print(p)
    
    if (length(sample_names) > 1) {
      cur_df$pval_start[is.na(cur_df$pval_start)] = 1
      cur_df[,'pval_start'] = as.numeric(as.character(cur_df[,'pval_start']))
      
      cur_df[((cur_df$pval_start <= SIGNIFICANCE_LEVELS[1])), 'pval'] = SIGNIFICANCE_LEVELS[1]
      cur_df[((cur_df$pval_start > SIGNIFICANCE_LEVELS[1]) & (cur_df$pval_start <= SIGNIFICANCE_LEVELS[2])), 'pval'] = SIGNIFICANCE_LEVELS[2]
      cur_df[((cur_df$pval_start > SIGNIFICANCE_LEVELS[2]) & (cur_df$pval_start <= SIGNIFICANCE_LEVELS[3])), 'pval'] = SIGNIFICANCE_LEVELS[3]
      cur_df[((cur_df$pval_start > SIGNIFICANCE_LEVELS[3]) & (cur_df$pval_start <= SIGNIFICANCE_LEVELS[4])), 'pval'] = SIGNIFICANCE_LEVELS[4]
      cur_df[cur_df$pval_start > SIGNIFICANCE_LEVELS[4], 'pval'] = SIGNIFICANCE_LEVELS[5]
      cur_df$pval = factor(cur_df$pval, levels=unique(sort(cur_df$pval, decreasing = T)))
      cur_df$sample = factor(cur_df$sample, levels=unique((cur_df$sample)))
      
      
      myBreaks = SIGNIFICANCE_LEVELS
      myColors = ZERO_TO_ONE_PVAL_SCALE
      
      cur_p = cur_df[cur_df[,1]!=sample_names[1], c('sample', 'bin', 'pval_start')]
      cur_p_wide = reshape(cur_p, idvar = "sample", timevar = 'bin', direction='wide')
      rows_n = as.character(cur_p_wide[,1])
      cur_p_wide = cur_p_wide[,2:ncol(cur_p_wide)]
      cur_p_wide = apply(cur_p_wide, 2, as.numeric)
      cur_p_wide = as.data.frame(cur_p_wide)
      
      # assign the colnames for the heatmap
      if (ncol(cur_p_wide) != 1) {
        cur_cols = as.character(seq(1,ncol(cur_p_wide),1))
      }
      else {
        cur_cols = c('1')
      }
      
      # if nothing is significant, do not make a heatmap
      if (!(setequal(unique(cur_p$pval_start), as.numeric(1)))) {
        p2 = pheatmap(cur_p_wide,
                      cluster_rows = F,
                      cluster_cols = F,
                      color = myColors, main = paste0(alphabet[i] ,": bonferroni corrected p-val (permissive scale)"), 
                      labels_col = cur_cols,
                      labels_row = rows_n, 
                      cellwidth = 7, cellheight = 7)
        
        print(p2)
        
        myColors = PVAL_COLORS
        
        p3 = pheatmap(cur_p_wide,
                      cluster_rows = F,
                      cluster_cols = F,
                      color = myColors, main = paste0(alphabet[i] ,": bonferroni corrected p-val (only significant)"), 
                      labels_col = as.character(seq(1,ncol(cur_p_wide),1)),
                      labels_row = rows_n, breaks = myBreaks,
                      cellwidth = 7, cellheight = 7)
        
        print(p3)
      }

      else {
        print('no value was significant (all=1), no output heatmap')
      }

      plot_list=list()
      plot_list[[1]] = p
      plot_list[[2]] = p2[[4]]
      plot_list[[3]] = p3[[4]]
      g = grid.arrange(arrangeGrob(grobs= plot_list,nrow=3,ncol=1))
      g

    }
  }
  dev.off()
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
}



parse_plot_colors = function(plot_colors, sample_names) {
  if (!is.null(plot_colors)) {
    plot_colors = unlist(strsplit(plot_colors, split = ","))
  }
  else {
    coul = brewer.pal(4, "PuOr") 
    plot_colors = colorRampPalette(coul)(length(sample_names))
  }
  return(plot_colors)
}



calculate_pvals_for_matching = function(lst_ordered_dfs_bins_as_columns, sample_names) {

  pvals_df = data.frame(matrix(NA, nrow = (length(sample_names)-1), ncol = (ncol(lst_ordered_dfs_bins_as_columns[[1]]))))
  rownames(pvals_df) = sample_names[2:length(sample_names)]
  colnames(pvals_df) = paste('bin', as.character(seq(1:ncol(lst_ordered_dfs_bins_as_columns[[1]]))), sep='')

  for (j in 2:length(lst_ordered_dfs_bins_as_columns)) {
    for (u in 1:ncol(lst_ordered_dfs_bins_as_columns[[1]])) {
      result = lapply(lst_ordered_dfs_bins_as_columns, "[", , u)
      pvals_df[sample_names[j], u] = wilcox.test(result[[1]], result[[j]], paired = T)$p.value
    }
  }
  return(pvals_df)
}



plot_repeats = function(df_ratios, df_pvals, sample_names, out_full_path) {
    
  df_pvals = apply(df_pvals, 1, function(x) {p.adjust(as.numeric(x), method = "bonferroni", n=ncol(df_pvals))})
  
  a = ggplot(df_ratios, aes(x=bin, y=rep_frequency, color=origin)) + 
    geom_line() + 
    geom_point() + 
    scale_color_discrete(labels = sample_names) + 
    theme_bw()
    
  b = pheatmap(t(df_pvals),
               cluster_rows = F,
               cluster_cols = F,
               labels_col = rownames(df_pvals),
               labels_row = colnames(df_pvals), 
               color = PVAL_COLORS, breaks = SIGNIFICANCE_LEVELS,
               main = paste0("bonferroni corrected p-val"),
               cellwidth = 10, cellheight = 10)

  pdf(out_full_path, width=14)
  print(a)
  grid::grid.newpage()
  print(b)
  dev.off()
}



calculate_lower_ratio = function(df, alphabet, sample_names, matching=NULL, plot_colors, out_full_path) {
  
  plot_colors = parse_plot_colors(plot_colors, sample_names)
  
  df = df[,c(1,3,10,13:ncol(df),11)]
  df = as.data.frame(apply(df, 2, as.character))
  
  alphabet = tolower(alphabet)
  
  lower_ratio = function(x) {sum(str_count(as.character(x), alphabet))/nchar(as.character(x))}
  # calc ratio:
  df[,c(3:ncol(df))] = as.data.frame(lapply(df[,c(3:ncol(df))], function(x) {sapply(x, lower_ratio)}))

  
  if (!is.null(matching)) {
    df$control_for = apply(df, 1, function(x) {ifelse(as.numeric(x['origin']) == 1, x['name'], matching[[as.numeric(x['origin'])-1]][matching[[as.numeric(x['origin'])-1]]$name == x['name'], 'control_for'])})
    split_var = list(df$origin, df$control_for)
    dfs = split(df, split_var)

    dfs_mean = lapply(dfs, function(x) {apply(x[3:(ncol(x)-1)], 2, mean)})
    
    combine_by_origin = list()
    for (j in 1:length(sample_names)) {
      combine_by_origin[[j]] = dfs_mean[startsWith(names(dfs_mean), as.character(j))]
      
    }
    dfs_back = lapply(combine_by_origin, function(x) {t(as.data.frame(x))})

    # taking into acocunt that the data frames are now ordered in the same matching order:
    mean_repeats = lapply(dfs_back, function(x) {apply(x, 2, mean)})
    mean_repeats_df = t(do.call(rbind, mean_repeats))
    df_lst = list()
    for (i in 1:length(sample_names)) {
      ratio_df = data.frame("origin" = c(rep(sample_names[i], nrow(mean_repeats_df))),
                            "rep_frequency" = as.numeric(mean_repeats_df[,i]))
      rownames(ratio_df) = NULL
      ratio_df$bin = as.numeric(rownames(ratio_df))
      df_lst[[i]] = ratio_df
    }
    
    df_comb = do.call(rbind, df_lst)
    df_comb$origin = factor(df_comb$origin)

    pvals_df = calculate_pvals_for_matching(dfs_back, sample_names)
  }
  
  else { # if (is.null(matching)) {
    dfs = split(df, df$origin)
    df_lst = list()

    for (i in 1:length(unique(df$origin))) {
      cur_df = dfs[[i]]
      cur_df_rep_ratio = apply(cur_df[,c(3:ncol(cur_df))], 2, mean)
      ratio_df = data.frame("origin" = c(rep(unique(df$origin)[i], length(cur_df_rep_ratio))),
                            "rep_frequency" = as.numeric(cur_df_rep_ratio))
      rownames(ratio_df) = NULL
      ratio_df$bin = as.numeric(rownames(ratio_df))
      df_lst[[i]] = ratio_df
    }

    df_comb = do.call(rbind, df_lst)
    df_comb$origin = factor(df_comb$origin)
    
    # calculate pvals:
    if (length(df_lst) > 1) {
      pvals_lst = list()

      for (m in 2:length(dfs)) {
        pvals = c()

        for (j in 3:ncol(dfs[[1]])) {
          cur_p = wilcox.test(dfs[[1]][,j], dfs[[m]][,j])$p.value
          pvals = c(pvals, cur_p)
        }

        pvals_lst[[m-1]] = pvals 
      }

    }
    pvals_df = do.call(rbind, pvals_lst)
       
  }

  plot_repeats(df_comb, pvals_df, sample_names, out_full_path)
  return(ratio_df)
}


create_alphabet_metagene = function(fasta_path_vec, match=NULL, alphabet_vec, sample_names_vec,
                                    output_title=NULL, out_full_path, make_upper, plot_colors) {
  
  alphabet_vec = parse_input_alphabet(alphabet_vec)
  fastas = read_fastas(fasta_path_vec, make_upper)
  
  if (is.null(sample_names_vec)) {
    sample_names_vec = paste0(rep('sample', length(fastas)), as.character(seq(1,length(fastas),1)))
  }
  else {
    sample_names_vec = unlist(strsplit(sample_names_vec, split = ","))
  }
  
  if (!is.null(match)) {
    matching = read_matching_control_data(match)
    if (!(is.list(matching))) {
      return('Wrong matching files. Exited program.')
    }
  }
  
  bound_fastas = bind_fastas(fastas)
  
  num_bins = define_num_bins(bound_fastas)
  combined_fasta = num_bins[[2]]
  num_bin = num_bins[[1]]
  
  comb_fasta = get_num_bases_per_bin(combined_fasta, num_bin)
  
  combined = divide_bins(comb_fasta)

  
  # next part is relevant if one wants to calculate the rate of repetitive elements:
  if (!(make_upper)) {
    calculate_lower_ratio(combined, alphabet_vec, sample_names_vec, matching, plot_colors, out_full_path)
    return('plotted lowercase alphabet overall ratio successfully- to plot metagene for alphabet remove -u option')
  }
  
  if (!is.null(match) & (length(sample_names_vec) > 1)) {
    wilcoxon_calc = calculate_wilcoxon(alphabet_vec, combined, sample_names_vec, matching, num_bin)
    bin_freq = calculate_bin_frequency(alphabet_vec, combined, matching)
    arranged_df = arrange_df_to_plot(bin_freq, alphabet_vec, sample_names_vec, wilcoxon_calc)
  }

  else if (is.null(match) & (length(sample_names_vec) > 1)) {
    wilcoxon_calc = calculate_wilcoxon(alphabet_vec, combined, sample_names_vec, num_bin = num_bin)
    bin_freq = calculate_bin_frequency(alphabet_vec, combined)
    arranged_df = arrange_df_to_plot(bin_freq, alphabet_vec, sample_names_vec, wilcoxon_calc)
  }

  else {
    bin_freq = calculate_bin_frequency(alphabet_vec, combined)
    arranged_df = arrange_df_to_plot(bin_freq, alphabet_vec, sample_names_vec, NULL)
  }
  
  melted = melt_df(arranged_df, sample_names_vec)
    
  print(melted)
  p = plot_metagene(melted, alphabet_vec, sample_names_vec, output_title, out_full_path, plot_colors)

  print("finished successfuly")
  print(paste0("saved pdf path is: ", out_full_path))
  
}



create_alphabet_metagene(opt$fasta_path_vec,
                         opt$match,
                         opt$alphabet_vec,
                         opt$sample_names_vec,
                         opt$output_title,
                         opt$out_full_path,
                         opt$make_upper,
                         opt$plot_colors)

warnings()

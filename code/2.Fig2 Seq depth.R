####################
#Set working directory
####################
path_to_wd = "path/to/code/dir" # path to the "code" folder
setwd(path_to_wd)

####################
#Load required libraries
####################
library(dplyr)
library(tidyverse)

library(ape)

library(ggplot2)
library(ggh4x) #facet_nested function and ggplot
library(cowplot)

source("miscellaneous funcs.R")

####################
#Set paths to file
####################
#inputs
variant_table_dat_filename <- "data/variant table.txt"
ref_annotations_filename <- "data/ref genome annotation.txt"

#outputs
formated_avg_seq_dep_dat_filename <- "results/Supplementary Table 1.avg_seq_dep_by_passage_stock.txt"
Fig2_seq_dep.png_filename <- "results/Fig2_seq_dep.png"
Fig2_seq_dep.svg_filename <- "results/Fig2_seq_dep.svg"

Misc_info_avg_seq_dep_filename <- "results/Misc_info.avg_seq_dep.txt"

####################
#Set parameter values
####################
depth_lw_threshold <- 50 

####################
#Make site-wise sequencing depth plots
####################
#load sequencing depth (and variant) data
variant_table_dat <- read.table(variant_table_dat_filename, header = T, sep = "\t") %>% 
	set_factors #from "miscellaneous funcs.R"

#compute average sequencing depths of each passage stock and write to file
variant_table_dat %>% 
	group_by(variant, sample, cell_line, replicate, passage) %>% 
	summarise(depth = round(mean(depth),2)) %>%
	spread(passage, depth) %>% as.data.frame %>%
	write.table(file = formated_avg_seq_dep_dat_filename, quote = F, sep = "\t", row.name = F)

#compute average sequencing depths and write to file
#across all datasets
avg_seq_dep_overall <- variant_table_dat %>% 
	summarise(depth = round(mean(depth),2)) %>% 
	ungroup()

#by sample
avg_seq_dep_by_sample <- variant_table_dat %>% 
	group_by(variant, sample) %>% 
	summarise(depth = round(mean(depth),2)) %>% 
	ungroup() %>% as.data.frame() 

sink(Misc_info_avg_seq_dep_filename)
cat("Average sequencing depths\n")
cat("=========================\n")
cat(sprintf("Across all datasets: %sx\n", avg_seq_dep_overall))
cat("By sample\n")
print(avg_seq_dep_by_sample)
cat("\n")
sink()

#compute average sequencing depths by position across all datasets
avg_dep_by_position <- variant_table_dat %>% 
	group_by(position) %>% 
	summarise(depth = round(mean(depth),2)) %>% 
	ungroup() %>% as.data.frame() %>%
	mutate(variant_sample = "All samples", .before = position)

#get info on genome regions with poor sequencing depths (<50x) and write it to file
ref_annotations <- read.table("data/ref genome annotation.txt", header = TRUE, sep = "\t", quote = "")

avg_dep_at_low_dep_position <- avg_dep_by_position %>% 
					filter(depth < depth_lw_threshold) #%>% filter(30 < position, position < 29675)

avg_dep_at_low_dep_position$gene <- get_gene_from_pos(
	avg_dep_at_low_dep_position$position,
	ref_annotations %>% 
		add_row(start = 1, end = 265, type = "UTR", gene = "5'UTR", .before = 1) %>% 
		add_row(start = 29675, end = 29903, type = "UTR", gene = "3'UTR")
)

sink(Misc_info_avg_seq_dep_filename, append = T)
cat(sprintf("Genome regions with poor sequencing depths (<%dx)\n", depth_lw_threshold))
cat("=========================\n")
avg_dep_at_low_dep_position %>% group_by(gene) %>% 
	summarise(
		range_position = paste(range(position), collapse = "-"),
		range_avg_dep = paste(range(depth), collapse = "-")
	) %>% ungroup() %>% as.data.frame()
sink()

#plot sequencing depths averaged across all datasets
seq_dep_plot_overall <- variant_table_dat %>% 
	mutate(
		variant_sample = "All samples", 
		depth = ifelse(depth == 0, 0.1, depth), 
		exp_id = paste(sample, cell_line, passage, replicate, sep = "_")
	) %>% 
	ggplot(aes(x = position, y = depth)) +
	geom_hline(yintercept = depth_lw_threshold, linetype = "dashed") + 
	geom_line(aes(group = exp_id), col = "grey", alpha = .5) + 
	geom_line(
		data = avg_dep_by_position %>% mutate(depth = ifelse(depth == 0, 0.1, depth)),
		col = "black",
		size = 0.5
	) +
	scale_x_continuous(name = "Position", limits = c(1,29903), expand = c(0, 0)) +
	scale_y_continuous(name = "Sequencing depth", trans = "log10", limits = c(1,NA), expand = c(0, 0))+
	facet_grid(variant_sample ~ ., 
		labeller = labeller(variant_sample = c("All samples" = "All\nsamples"))
	) +
	theme_classic() +
	labs(title = "Overall") +
	theme(
		legend.position = "none",
		plot.title = element_text(face = "bold", size = 10),
		axis.title.x = element_blank(),
		axis.text.x = element_blank(), #axis.ticks.x = element_blank(),
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank()
	)

#compute average sequencing depths by virus
avg_dep_by_virus_and_cell_line <- variant_table_dat %>% 
	group_by(variant, sample, variant_sample, cell_line, position) %>% 
	summarise(depth = mean(depth)) %>%
	ungroup() %>% as.data.frame()

#plot sequencing depths averaged by virus
seq_dep_plot_by_virus_and_cell_line <- variant_table_dat %>% 
	mutate(
		depth = ifelse(depth == 0, 0.1, depth),
		exp_id = paste(sample, cell_line, passage, replicate, sep = "_")
	) %>% 
	ggplot(aes(x = position, y = depth)) +
	geom_hline(yintercept = depth_lw_threshold, linetype = "dashed") + 
	geom_line(aes(group = exp_id), col = "grey", alpha = .5) + 
	geom_line(
		data = avg_dep_by_virus_and_cell_line %>% mutate(depth = ifelse(depth == 0, 0.1, depth)),
		aes(col = variant_sample, linetype = variant_sample), 
		size = 0.5
	) +

	scale_color_manual(
		name = "Virus sample",
		values = c(
			"B.1.36.16: 73NLt" = "#00A087FF", 
			"B.1.36.16: CV130" = "#91D1C2DD", 
			"AY.30: OTV54" = "#E64B35FF", 
			"AY.30: NH783" = "#F39B7FFF" 
		)
	) +
	scale_linetype_manual(
		name = "Virus sample",
		values = c( 
			"B.1.36.16: 73NLt" = "solid", 
			"B.1.36.16: CV130" = "solid", 
			"AY.30: OTV54" = "solid", 
			"AY.30: NH783" = "solid" 
		)
	) +
	
	scale_x_continuous(name ="Position", limits = c(1,29903), expand = c(0, 0)) +
	scale_y_continuous(name ="Sequencing depth", trans = "log10", limits = c(1,NA), expand = c(0, 0)) +
	facet_nested(
		variant + sample + cell_line ~ . ,
		labeller = labeller(cell_line = c("Clinical sample" = "Clinical\nsample", "Vero E6" = "Vero E6", "Vero E6-TMPRSS2" = "Vero E6/\nTMPRSS2", "Calu-3" = "Calu-3"))
	) +

	theme_classic() +
	labs(title = "By virus and cell line") +
	theme(
		legend.position = "none",
		plot.title = element_text(face = "bold", size = 10),
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank()
	)

####################
#Make a composit plot of genome structure, and sequencing depth
####################
aligned_left_plots <- align_plots(
	make_ref_genome_structure_plot(ref_annotations_filename = ref_annotations_filename),
	seq_dep_plot_overall,
	seq_dep_plot_by_virus_and_cell_line,
	align = "v", axis = "l"
)

Fig2_plot <- align_plots(
	aligned_left_plots[[1]],
	aligned_left_plots[[2]],
	aligned_left_plots[[3]],
	align = "v", axis = "lr"
)

Fig2_plot <- ggdraw() + 
	draw_plot(Fig2_plot[[3]], x = 0, y = 0.00, width = 1.0, height = 0.88) +
	draw_plot(Fig2_plot[[2]], x = 0, y = 0.88, width = 1.0, height = 0.08) +
	draw_plot(Fig2_plot[[1]], x = 0, y = 0.95, width = 1.0, height = 0.05)

Fig2_plot <- Fig2_plot +
	theme(plot.background = element_rect(fill = "white", color = NA),  panel.border = element_blank())

ggsave(Fig2_seq_dep.png_filename,
	plot = Fig2_plot,
	width = 20, height = 32, units = "cm",
	dpi = 300)

ggsave(Fig2_seq_dep.svg_filename,
	plot = Fig2_plot,
	width = 20, height = 32, units = "cm",
	dpi = 300)


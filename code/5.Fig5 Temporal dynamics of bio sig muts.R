####################
#Set working directory
####################
path_to_wd = "path/to/code/dir" # path to the "code" folder
setwd(path_to_wd)

####################
#Load required libraries
####################
library(dplyr)

library(ggplot2)
library(ggh4x) #facet_nested function
library(cowplot) #align_plots

####################
#Set paths to file
####################
#inputs
variant_table_dat_MBCS_filename <- "data/changes at MBCS.txt"
variant_table_dat_SH665Y_filename <- "data/changes at S655.txt"
variant_table_dat_NLinker_filename <- "data/changes at N linker.txt"
variant_table_dat_NSP1_filename <- "data/changes at NSP1.txt"
variant_table_dat_M_filename <- "data/changes at M.txt"
variant_table_dat_NSP3_PLpro_filename <- "data/changes at NSP3 PLpro.txt"

#outputs
Fig5_Tem_dyn_of_bio_sig_muts.raw.png_filename <- "results/Fig5_Temporal dynamics of bio sig muts.raw.png"
Fig5_Tem_dyn_of_bio_sig_muts.raw.svg_filename <- "results/Fig5_Temporal dynamics of bio sig muts.raw.svg"

####################
#Make individual mutation frequency over time plots
####################
#function for plotting
plot_temp_dyn_of_adap_change <- function(variant_table_dat, mut_pattern_levels){
	variant_table_dat <- variant_table_dat %>% 
	mutate(
		variant = factor(variant, levels = c("B.1.36.16", "AY.30")),
		sample = factor(sample, levels = c("73NLt", "CV130", "OTV54", "NH783")),
		cell_line = factor(cell_line, c("Vero E6", "Vero E6/TMPRSS2", "Calu-3")),
		passage = factor(passage, levels = c("P0", "P1", "P2", "P3", "P4"), ordered = TRUE),
		replicate = paste("Replicate: ", replicate, sep = ""),
		replicate = factor(replicate, levels = c("Replicate: A", "Replicate: B")),
		mut_pattern = factor(mut_pattern, levels = mut_pattern_levels), #the first one of the mut_pattern_levels is the WT one
		frequency = as.numeric(frequency),

		variant_sample = factor(paste(variant, sample, sep = ": "), levels = c("B.1.36.16: 73NLt", "B.1.36.16: CV130", "AY.30: OTV54", "AY.30: NH783")) #add variant_sample column
	)

	cols_dict <- setNames(
		c("#FFFFFFFF", rainbow(variant_table_dat %>% select(mut_pattern) %>% n_distinct() - 1, alpha = 1)),
		variant_table_dat %>% 
		group_by(mut_pattern) %>%
		summarise(sum_freq = sum(frequency)) %>% pull(mut_pattern)
	); cols_dict["Others"] <- "#808080FF"

	mu_freq_plots <- 
		ggplot(data = variant_table_dat, aes(fill = mut_pattern, y = frequency, x = passage)) + 
		geom_bar(position = "fill", stat = "identity", colour = "black", size = 0.1) +
		scale_fill_manual(name = "Mutational\npatterns", values = cols_dict) +

		scale_x_discrete(name = NULL) + scale_y_continuous(name = "Mutation\nfrequency") +
		facet_nested(cell_line ~ position + variant_sample + replicate) +

		theme_classic() +
		theme(
			plot.title = element_text(size = 10, face = "bold"),
			plot.margin = unit(c(0,2.5,0,0),"cm"),

			legend.justification = c(0, 0),
			legend.position = c(1.05, 0.0),
			legend.title = element_text(size = 8, face = "bold"), 
			legend.text = element_text(size = 6, family = "mono"),
			legend.key.width = unit(0.25,"cm"), 
			legend.key.height = unit(0.25,"cm"), 
			legend.box.margin = margin(0, 0, 0, 0),
			legend.margin = margin(t = 0, b = 0, unit = "cm"),

			strip.background = element_rect(fill = "white", color = "black"),
			strip.text = element_text(size = 6, colour = "black"),

			axis.text = element_text(size = 6),
			axis.title = element_text(size = 8, face = "bold"),
		) +
		guides(fill = guide_legend(ncol = 1))

	return(mu_freq_plots)
}

#load variant datasets
variant_table_dat_MBCS <- read.table(variant_table_dat_MBCS_filename, header = T, sep = "\t", quote = "")
variant_table_dat_MBCS_B.1.36.16 <- variant_table_dat_MBCS %>% filter(variant == "B.1.36.16") 
variant_table_dat_MBCS_AY.30 <- variant_table_dat_MBCS %>% filter(variant == "AY.30")

variant_table_dat_SH665Y <- read.table(variant_table_dat_SH665Y_filename, header = T, sep = "\t", quote = "")

variant_table_dat_NLinker <- read.table(variant_table_dat_NLinker_filename, header = T, sep = "\t", quote = "")
variant_table_dat_NLinker_B.1.36.16 <- variant_table_dat_NLinker %>% filter(variant == "B.1.36.16") 
variant_table_dat_NLinker_AY.30 <- variant_table_dat_NLinker %>% filter(variant == "AY.30")

variant_table_dat_NSP1 <- read.table(variant_table_dat_NSP1_filename, header = T, sep = "\t", quote = "")

variant_table_dat_M <- read.table(variant_table_dat_M_filename, header = T, sep = "\t", quote = "")

variant_table_dat_NSP3_PLpro <- read.table(variant_table_dat_NSP3_PLpro_filename, header = T, sep = "\t", quote = "")

#plots individual dynamics
mu_freq_plots_MBCS_B.1.36.16 <- plot_temp_dyn_of_adap_change(
	variant_table_dat = variant_table_dat_MBCS_B.1.36.16, 
	mut_pattern_levels = c("YQTQTNSPRRAR", "Y-----SPRRAR", "YQTQTNSPWRAR", "YQTQTNSPQRAR", "YQTQTNSPLRAR", "Others")
	) +
	theme(
		plot.margin = unit(c(0,2.25,0,0),"cm"),
	)

mu_freq_plots_MBCS_AY.30 <- plot_temp_dyn_of_adap_change(
	variant_table_dat = variant_table_dat_MBCS_AY.30, 
	mut_pattern_levels = c("TNSRRRARSVA", "TNSRQRARSVA", "T-------SVA", "T----------", "TNSRWRARSVA", "TNSRRRARIVA", "TNSRRRARRVA", "TNSRRRACSVA", "TR---RARSVA", "TNSRRRARSVT", "TNSRQQARSVA", "Others")
	) +
	theme(
		plot.margin = unit(c(0,2.25,0,0),"cm"),
	)

mu_freq_plots_SH655Y <- plot_temp_dyn_of_adap_change(
	variant_table_dat = variant_table_dat_SH665Y, 
	mut_pattern_levels = c("H655", "Y655", "Others")
	) +
	theme(
		plot.margin = unit(c(0,1.5, 0,0),"cm"),
		#legend.position = c(1.075, 0.0),
	)

mu_freq_plots_NLinker_B.1.36.16 <- plot_temp_dyn_of_adap_change(
	variant_table_dat = variant_table_dat_NLinker_B.1.36.16, 
	mut_pattern_levels = c("A208", "V208", "Others")
	) +
	scale_y_continuous(name = NULL) +
	theme(
		plot.margin = unit(c(0,1.5,0,0),"cm"),
		#legend.position = c(1.30, 0.0)
	)
mu_freq_plots_NLinker_AY.30 <- plot_temp_dyn_of_adap_change(
	variant_table_dat = variant_table_dat_NLinker_AY.30,
	mut_pattern_levels = c("SMGTSPA", "S------", "S----->", "SMGPSPA", "Others")
	) +
	scale_y_continuous(name = NULL) +
	theme(
		plot.margin = unit(c(0,1.5,0,0),"cm"),
		#legend.position = c(1.30, 0.0)
	)

mu_freq_plots_NSP1 <- plot_temp_dyn_of_adap_change(
	variant_table_dat = variant_table_dat_NSP1,
	mut_pattern_levels = c("GHVMV", "---->", "GH--V", "GHV-V", "V---V", "GHV->", "Others")
	) +
	theme(
		plot.margin = unit(c(0,2.25,0,0),"cm"),
	)

mu_freq_plots_M <- plot_temp_dyn_of_adap_change(
	variant_table_dat = variant_table_dat_M,
	mut_pattern_levels = c("CLVGLMWLSYFIA", "------------>", "-------------", "CLV----------", "Others")
	) +
	theme(
		plot.margin = unit(c(0,2.25,0,0),"cm"),
	)

mu_freq_plots_NSP3_PLpro <- plot_temp_dyn_of_adap_change(
	variant_table_dat = variant_table_dat_NSP3_PLpro,
	mut_pattern_levels = c("S848", "C848", "Others")
	) +
	scale_x_discrete(name = "Passage") +
	theme(
		plot.margin = unit(c(0,2.25,0,0),"cm"),
	)

####################
#Make a raw composit plot, and write to file
####################
left_aligned_plots <- align_plots(
mu_freq_plots_MBCS_B.1.36.16,
mu_freq_plots_MBCS_AY.30,
mu_freq_plots_SH655Y,
mu_freq_plots_NSP1,
mu_freq_plots_M,
mu_freq_plots_NSP3_PLpro,
	align = "v", axis = "l")

mu_freq_plots_MBCS_B.1.36.16 <- left_aligned_plots[[1]]
mu_freq_plots_MBCS_AY.30 <- left_aligned_plots[[2]]
mu_freq_plots_SH655Y <- left_aligned_plots[[3]]
mu_freq_plots_NSP1 <- left_aligned_plots[[4]]
mu_freq_plots_M <- left_aligned_plots[[5]]
mu_freq_plots_NSP3_PLpro <- left_aligned_plots[[6]]

right_aligned_plots <- align_plots(
mu_freq_plots_MBCS_B.1.36.16,
mu_freq_plots_MBCS_AY.30,
mu_freq_plots_NLinker_AY.30,
mu_freq_plots_NSP1,
mu_freq_plots_M,
mu_freq_plots_NSP3_PLpro,
	align = "v", axis = "r")

mu_freq_plots_MBCS_B.1.36.16 <- right_aligned_plots[[1]]
mu_freq_plots_MBCS_AY.30 <- right_aligned_plots[[2]]
mu_freq_plots_NLinker_AY.30 <- right_aligned_plots[[3]]
mu_freq_plots_NSP1 <- right_aligned_plots[[4]]
mu_freq_plots_M <- right_aligned_plots[[5]]
mu_freq_plots_NSP3_PLpro <- right_aligned_plots[[6]]

mu_freq_plots_SH655Y_and_NLinker <- plot_grid(
	mu_freq_plots_SH655Y,
	mu_freq_plots_NLinker_B.1.36.16, 
	mu_freq_plots_NLinker_AY.30,
	labels = c("c)", "d)", "e)"),
	nrow = 1, rel_widths = c(1.55, 1, 1.2)
) 

mu_freq_plots <- plot_grid(
mu_freq_plots_MBCS_B.1.36.16,
mu_freq_plots_MBCS_AY.30,
mu_freq_plots_SH655Y_and_NLinker,
mu_freq_plots_NSP1,
mu_freq_plots_M,
mu_freq_plots_NSP3_PLpro,
	labels = c("a)", "b)", "", "f)", "g)", "h)"),
	nrow = 6, rel_heights = c(1, 1, 1, 1, 1, 2.0338)
) 

#save figure to files
ggsave(Fig5_Tem_dyn_of_bio_sig_muts.raw.png_filename, #need a bit of manual curation to make it pretty
	plot = mu_freq_plots,
	width = 16, height = 24, units = "cm",
	dpi = 300)

ggsave(Fig5_Tem_dyn_of_bio_sig_muts.raw.svg_filename, #need a bit of manual curation to make it pretty
	plot = mu_freq_plots,
	width = 16, height = 24, units = "cm",
	dpi = 300)

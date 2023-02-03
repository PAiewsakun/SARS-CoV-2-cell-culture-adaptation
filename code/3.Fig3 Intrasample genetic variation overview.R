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

library(ggplot2)
library(ggh4x) #facet_nested function and ggplot
library(cowplot)

source("miscellaneous funcs.R")

####################
#Set paths to file
####################
#inputs
variant_table_dat_filename <- "data/variant table.txt"

#outputs
Misc_info_var_site_count_by_cellline_filename <- "results/Misc_info.var_site_count_by_cell_line.txt"
Misc_info_var_site_count_dynamics_filename <- "results/Misc_info.var_site_count_dynamics.txt"

formated_const_and_var_site_counts_filename <- "results/Supplementary Table 2.const_and_var_site_counts.txt"
var_site_count_by_gen_change_type_filename <- "results/Supplementary Table 3.variant site count by genetic change type.txt"

Fig3_genetic_change_count.png_filename <- "results/Fig3_genetic_change_count.png"
Fig3_genetic_change_count.svg_filename <- "results/Fig3_genetic_change_count.svg"

####################
#Set parameter values
####################
mu_freq_lw_threshold <- 0.01

####################
#Identify constant sites, and for variant sites, count the number of variant types, averaged across passage stocks and experimental replicates
####################
#load variant data
variant_table_dat <- read.table(variant_table_dat_filename, header = T, sep = "\t") %>% 
	set_factors %>% #from "miscellaneous funcs.R"
	filter(30 < position, position < 29675) %>%
	mutate(
		fA = A/depth,
		fT = T/depth,
		fC = C/depth,
		fG = G/depth,
		fINS = INS/depth,
		fDEL = DEL/depth
	) %>% replace(is.na(.), 0)

#compute the number of variant types for each site, mutation freq cutouff == mu_freq_lw_threshold
variant_type_count_dat <- variant_table_dat %>%
	gather(gen.variant, freq, fA, fT, fC, fG, fINS, fDEL) %>% 
	group_by(variant, sample, cell_line, passage, replicate, position) %>%
	summarise(variant_types = sum(freq > mu_freq_lw_threshold)) %>% 
	ungroup() %>% as.data.frame() 

#count constant and variant sites
var_type_count_table <- variant_type_count_dat %>% 
	group_by(variant, sample, cell_line, passage, replicate) %>%
	summarise(
		const_sites = sum(variant_types == 1),
		var_sites = sum(variant_types > 1),
		bivar_sites = sum(variant_types == 2),
		multivar_sites = sum(variant_types > 2),
	) %>% 
	arrange(cell_line, variant, sample, passage, replicate) %>% 
	ungroup() %>% as.data.frame() 

#format the var_type_count_table and write it to file
var_type_count_table %>%
mutate(
	total_sites = const_sites + bivar_sites + multivar_sites,
	fconst_sites = round(const_sites/total_sites*100, 2),
	fvar_sites = round(var_sites/total_sites*100, 2),
	fbivar_sites = round(bivar_sites/var_sites*100, 2),
	fmultivar_sites = round(multivar_sites/var_sites*100, 2),

	const_sites = paste(prettyNum(const_sites, big.mark = ","), paste("(", fconst_sites, "%)", sep = ""), sep = " "),
	var_sites = paste(prettyNum(var_sites, big.mark = ","), paste("(", fvar_sites, "%)", sep = ""), sep = " "),
	bivar_sites = paste(prettyNum(bivar_sites, big.mark = ","), paste("(", fbivar_sites, "%)", sep = ""), sep = " "),
	multivar_sites = paste(prettyNum(multivar_sites, big.mark = ","), paste("(", fmultivar_sites, "%)", sep = ""), sep = " "),
) %>% 
select(cell_line, variant, sample, everything(), -(total_sites:fmultivar_sites)) %>%
write.table(formated_const_and_var_site_counts_filename, row.name = F, quote = F, sep = "\t")

#compute simple stats regarding variant site numbers and write to file
variant_site_count_dat <- variant_type_count_dat %>% 
	filter(variant_types >= 2) %>% 
	group_by(cell_line, variant, sample, passage, replicate) %>%
	summarise(variant_site_count = n()) %>% as.data.frame

#Clinical samples
variant_site_count_dat_clinical_samples <- variant_site_count_dat %>% 
	filter(cell_line == "Clinical sample") %>%
	group_by(cell_line) %>%
	summarise(
		med = median(variant_site_count),
		sd = round(sd(variant_site_count),2),
		range = paste(min(variant_site_count), "-", max(variant_site_count), sep = ""),
	)%>% as.data.frame

#Cultured samples
variant_site_count_dat_cultured_samples <- variant_site_count_dat %>% 
	filter(cell_line != "Clinical sample") %>% 
	mutate(cell_line = "Cultured sample") %>%
	group_by(cell_line) %>%
	summarise(
		med = median(variant_site_count),
		sd = round(sd(variant_site_count),2),
		range = paste(min(variant_site_count), "-", max(variant_site_count), sep = ""),
	)%>% as.data.frame

#by cell lines
variant_site_count_dat_by_call_line <- variant_site_count_dat %>% 
	filter(cell_line != "Clinical sample") %>%
	group_by(cell_line) %>%
	summarise(
		med = median(variant_site_count),
		sd = round(sd(variant_site_count),2),
		range = paste(min(variant_site_count), "-", max(variant_site_count), sep = ""),
	)%>% as.data.frame

#% of variant sites being bivariant in original clinical samples
fbivar_sites_clinical_samples <- var_type_count_table %>% 
	filter(cell_line == "Clinical sample") %>%
	mutate(fbivar_sites = bivar_sites/var_sites*100) %>%
	summarise(
		med = round(median(fbivar_sites), 2),
		sd = round(sd(fbivar_sites), 2),
		range = paste(round(min(fbivar_sites), 2), "-", round(max(fbivar_sites), 2), sep = "")
	)

#% of variant sites being bivariant in cultured samples
fbivar_sites_cultured_samples <- var_type_count_table %>%
	filter(cell_line != "Clinical sample") %>%
	mutate(fbivar_sites = bivar_sites/var_sites*100) %>%
	summarise(
		med = round(median(fbivar_sites), 2),
		sd = round(sd(fbivar_sites), 2),
		range = paste(round(min(fbivar_sites), 2), "-", round(max(fbivar_sites), 2), sep = "")
	)

#write the stats to file
sink(Misc_info_var_site_count_by_cellline_filename)
cat("Variant site numbers\n")
cat("=========================\n")
cat("Overall\n")
cat(sprintf("%s: median = %s sites; sd = %s sites; range = %s sites\n", variant_site_count_dat_clinical_samples$cell_line, variant_site_count_dat_clinical_samples$med, variant_site_count_dat_clinical_samples$sd, variant_site_count_dat_clinical_samples$range))
cat(sprintf("%s: median = %s sites; sd = %s sites; range = %s sites\n", variant_site_count_dat_cultured_samples$cell_line, variant_site_count_dat_cultured_samples$med, variant_site_count_dat_cultured_samples$sd, variant_site_count_dat_cultured_samples$range))
cat("\n")
cat("By cell line\n")
for(i in 1:nrow(variant_site_count_dat_by_call_line)){
	cat(sprintf("%s: median = %s sites; sd = %s sites; range = %s sites\n", variant_site_count_dat_by_call_line[i, "cell_line"], variant_site_count_dat_by_call_line[i, "med"], variant_site_count_dat_by_call_line[i, "sd"], variant_site_count_dat_by_call_line[i, "range"]))
}
cat("\n")
cat("Proportion of variant sites being bivariant\n")
cat("=========================\n")
cat(sprintf("Clinical sample: median %s%%; range %s%%\n", fbivar_sites_clinical_samples$med, fbivar_sites_clinical_samples$range))
cat(sprintf("Cultured sample: median %s%%; range %s%%\n", fbivar_sites_cultured_samples$med, fbivar_sites_cultured_samples$range))
cat("\n")
sink()

#fit model to the polymoprhic site count to see how it varies with time
x <- var_type_count_table %>% filter(cell_line != "Clinical sample") %>% 
mutate(
	passage = ifelse(passage == "P1", 0, ifelse(passage == "P2", 1, ifelse(passage == "P3", 2, ifelse(passage == "P4", 3, NA)))), 
	variant_sample = paste(variant, sample, sep = "_")
	) 

mdl_0 <- lm(data = x, formula = var_sites ~ cell_line + passage)
mdl_1 <- lm(data = x, formula = var_sites ~ cell_line * passage)
anova(mdl_0, mdl_1)

mdl_2 <- lm(data = x, formula = var_sites ~ cell_line + passage + sample)
anova(mdl_0, mdl_2)

#write the stats to file
sink(Misc_info_var_site_count_dynamics_filename)
cat("Dynamics of variant site numbers through passages\n")
cat("=========================\n")
cat("H0: var_sites ~ cell_line + passage\n")
summary(mdl_0)
cat("H1: var_sites ~ cell_line * passage\n")
summary(mdl_1)
cat("H2: var_sites ~ cell_line + passage + sample\n")
summary(mdl_2)

cat("The rates at which the variant site count decreases are not significantly different among cell lines used in the propagation\n")
cat("=========================\n")
anova(mdl_0, mdl_1)
cat("\n")
cat("The rates at which the variant site count decreases are not significantly different among viral samples\n")
cat("=========================\n")
anova(mdl_0, mdl_2)
cat("\n")
sink()

#####################
#Plot polymorphic site counts through passages 
#####################
var_count_through_passages_plot <- 
	ggplot() + 
	geom_hline(
		data = var_type_count_table %>% 
			filter(cell_line == "Clinical sample") %>% 
			mutate(
				variant_sample = factor(paste(variant, sample, sep = ": "), levels = c("B.1.36.16: 73NLt", "B.1.36.16: CV130", "AY.30: OTV54", "AY.30: NH783"))
			) %>%
			select (-cell_line), 
		aes(yintercept = var_sites, col = variant_sample), linetype = "solid", lwd = 0.5) + 

	geom_point(
		data = var_type_count_table %>% filter(cell_line != "Clinical sample") %>%
		mutate(
			variant_sample = factor(paste(variant, sample, sep = ": "), levels = c("B.1.36.16: 73NLt", "B.1.36.16: CV130", "AY.30: OTV54", "AY.30: NH783")),
			exp_id = paste(sample, cell_line, replicate, sep = "_")
		), 
		aes(x = passage, y = var_sites, col = variant_sample), size = 1) + 

	geom_line(
		data = var_type_count_table %>% filter(cell_line != "Clinical sample") %>%
		mutate(
			variant_sample = factor(paste(variant, sample, sep = ": "), levels = c("B.1.36.16: 73NLt", "B.1.36.16: CV130", "AY.30: OTV54", "AY.30: NH783")),
			exp_id = paste(sample, cell_line, replicate, sep = "_")
		),
		aes(x = passage, y = var_sites, group = exp_id, col = variant_sample, linetype = replicate), linewidth = 0.5) +

	scale_color_manual(
		name = "Virus sample",
		values = c(
			"B.1.36.16: 73NLt" = "#00A087FF", 
			"B.1.36.16: CV130" = "#91D1C2DD", 
			"AY.30: OTV54" = "#E64B35FF", 
			"AY.30: NH783" = "#F39B7FFF" 
		),
		guide = guide_legend(title.position = "top", direction = "horizontal", keyheight = 0.4, nrow = 2, order = 1)

	) +
	scale_linetype_manual(
		name = "Replicate",
		values = c( 
			"A" = "solid", 
			"B" = "dashed"
		), 
		guide = guide_legend(title.position = "top", direction = "horizontal", keyheight = 0.4, order = 2)
	) +

	scale_x_discrete(name = "Passage") +
	scale_y_continuous(name = "Polymorphic site count") +

	geom_line(
		data = mdl_0$model %>% 
				select(cell_line, passage) %>%
				unique %>%
				mutate(var_sites = predict(mdl_0, data.frame(cell_line = cell_line, passage = passage))) %>%
				mutate(passage = recode(passage, "0" = "P1", "1" = "P2", "2" = "P3", "3" = "P4")) %>%
				mutate(exp_id = cell_line),
		aes(x = passage, y = var_sites, group = exp_id), col = "black", linetype = "solid", linewidth = 1) + 

	facet_grid(
		. ~ cell_line, 
		labeller = labeller(cell_line = c("Vero E6" = "Vero E6", "Vero E6-TMPRSS2" = "Vero E6/\nTMPRSS2", "Calu-3" = "Calu-3"))
	) +

	theme_light() + 
	theme(
		legend.justification = c(0, 1),
		legend.position = c(0.01, 0.98),
		#legend.position = "bottom",
		legend.background = element_rect(fill = alpha("white", 0.5)),
		legend.key = element_rect(fill = alpha("white", 0.5)),
		legend.box = "vertical",

		legend.title = element_text(size = 8, face = "bold"), 
		legend.text = element_text(size = 6),

		legend.margin = margin(0,0,0,0),
		legend.box.margin = margin(0,0,0,0),
		legend.spacing.y = unit(0.0, 'cm'), 

		axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6),

		strip.background =element_rect(fill = "white", color = "black"),
		strip.text = element_text(size = 8, color = "black", face = "bold")
	)

#####################
#Plot distributions of polymorphic site counts by genetic change type
#####################
#get site-wise original variants of the virus samples - major variants in the original samples
ori.var <- variant_table_dat %>% 
	filter(cell_line == "Clinical sample") %>% 
	mutate(ori.var = c("A", "T", "C", "G", "INS", "DEL")[max.col(select(., A, T, C, G, INS, DEL), ties.method = "first")]) %>%
	mutate(ori.var = ifelse(depth == 0, as.character(REF), ori.var)) %>% #if depth == 0, use REF var. REF is a factor... need as.character, otherwise will be number
	select(sample, position, ori.var) %>%
	arrange(sample, position)

#add ori.var to variant_table_dat
variant_table_dat <- left_join(variant_table_dat, ori.var, by = c("sample" = "sample", "position" = "position")) %>% 
	arrange(variant, sample, cell_line, position) 

#compute original variant frequency 
variant_table_dat$fOri.var <- variant_table_dat[1:nrow(variant_table_dat), c("fA", "fT", "fC", "fG", "fINS", "fDEL")][cbind(1:nrow(variant_table_dat), paste("f", variant_table_dat$ori.var, sep = ""))]

#determine genetic changes' frequencies and the number of types
x <- variant_table_dat %>% select(fA:fDEL)
x[x < mu_freq_lw_threshold] <- 0
variant_table_dat$fMut <- rowSums(x) - variant_table_dat$fOri.var

x <- variant_table_dat %>% select(fA:fDEL)
variant_table_dat$var.type_count <- rowSums(x >= mu_freq_lw_threshold) - 1 
variant_table_dat <- variant_table_dat %>% 
	mutate(var.type_count = ifelse(var.type_count < 0, 0, var.type_count))

#determine genetic changes' direction. For multivariant sites, the direction is weighted by alternative variant freq.
variant_table_dat <- variant_table_dat %>% 
	mutate(
		A.T = ifelse(ori.var == "A" & fT > mu_freq_lw_threshold, fT/fMut, 0),
		A.C = ifelse(ori.var == "A" & fC > mu_freq_lw_threshold, fC/fMut, 0),
		A.G = ifelse(ori.var == "A" & fG > mu_freq_lw_threshold, fG/fMut, 0),

		T.A = ifelse(ori.var == "T" & fA > mu_freq_lw_threshold, fA/fMut, 0),
		T.C = ifelse(ori.var == "T" & fC > mu_freq_lw_threshold, fC/fMut, 0),
		T.G = ifelse(ori.var == "T" & fG > mu_freq_lw_threshold, fG/fMut, 0),

		C.A = ifelse(ori.var == "C" & fA > mu_freq_lw_threshold, fA/fMut, 0),
		C.T = ifelse(ori.var == "C" & fT > mu_freq_lw_threshold, fT/fMut, 0),
		C.G = ifelse(ori.var == "C" & fG > mu_freq_lw_threshold, fG/fMut, 0),

		G.A = ifelse(ori.var == "G" & fA > mu_freq_lw_threshold, fA/fMut, 0),
		G.T = ifelse(ori.var == "G" & fT > mu_freq_lw_threshold, fT/fMut, 0),
		G.C = ifelse(ori.var == "G" & fC > mu_freq_lw_threshold, fC/fMut, 0),

		A.INS = ifelse(ori.var == "A" & fINS > mu_freq_lw_threshold, fINS/fMut, 0),
		T.INS = ifelse(ori.var == "T" & fINS > mu_freq_lw_threshold, fINS/fMut, 0),
		C.INS = ifelse(ori.var == "C" & fINS > mu_freq_lw_threshold, fINS/fMut, 0),
		G.INS = ifelse(ori.var == "G" & fINS > mu_freq_lw_threshold, fINS/fMut, 0),

		A.DEL = ifelse(ori.var == "A" & fDEL > mu_freq_lw_threshold, fDEL/fMut, 0),
		T.DEL = ifelse(ori.var == "T" & fDEL > mu_freq_lw_threshold, fDEL/fMut, 0),
		C.DEL = ifelse(ori.var == "C" & fDEL > mu_freq_lw_threshold, fDEL/fMut, 0),
		G.DEL = ifelse(ori.var == "G" & fDEL > mu_freq_lw_threshold, fDEL/fMut, 0),

		INS.A = ifelse(ori.var == "INS" & fA > mu_freq_lw_threshold, fA/fMut, 0),
		INS.T = ifelse(ori.var == "INS" & fT > mu_freq_lw_threshold, fT/fMut, 0),
		INS.C = ifelse(ori.var == "INS" & fC > mu_freq_lw_threshold, fC/fMut, 0),
		INS.G = ifelse(ori.var == "INS" & fG > mu_freq_lw_threshold, fG/fMut, 0),

		DEL.A = ifelse(ori.var == "DEL" & fA > mu_freq_lw_threshold, fA/fMut, 0),
		DEL.T = ifelse(ori.var == "DEL" & fT > mu_freq_lw_threshold, fT/fMut, 0),
		DEL.C = ifelse(ori.var == "DEL" & fC > mu_freq_lw_threshold, fC/fMut, 0),
		DEL.G = ifelse(ori.var == "DEL" & fG > mu_freq_lw_threshold, fG/fMut, 0),

		INS.DEL = ifelse(ori.var == "INS" & fDEL > mu_freq_lw_threshold, fDEL/fMut, 0),
		DEL.INS = ifelse(ori.var == "DEL" & fINS > mu_freq_lw_threshold, fINS/fMut, 0)
	)

#count genetic changes
genetic_change_count_table <- variant_table_dat %>% 
	group_by(variant, sample, cell_line, passage, replicate) %>%
	summarise(across(A.T:DEL.INS, sum)) %>% 
	arrange(cell_line) %>% 
	ungroup() %>% as.data.frame()

#write change counts by genetic type to file
genetic_change_count_table %>% write.table(var_site_count_by_gen_change_type_filename, row.name = F, quote = F, sep = "\t")

#plot variant site counts
genetic_change_count_table_long <- genetic_change_count_table %>% 
	gather(key = "gen.change", value = "count", A.T:DEL.G) %>% 
	mutate(exp_id = paste(variant, sample, sep = ": ")) 

genetic_change_count_table_long_mean_sd <- genetic_change_count_table_long %>%
	group_by(variant, sample, cell_line, gen.change) %>%
	summarise(
		se = sd(count, na.rm = TRUE)/sqrt(sum(!is.na(count))),
		count = mean(count),
		count_lb = count - se,
		count_ub = count + se,
	) %>% ungroup() %>% as.data.frame %>%
	mutate(exp_id = paste(variant, sample, sep = ": ")) 

genetic_change_count_plot <- genetic_change_count_table_long_mean_sd %>%
	ggplot(aes(x = gen.change, y = count, fill = exp_id, color = exp_id, linetype = exp_id)) +
	geom_bar(stat = "identity", width = 0.75) +
	geom_errorbar(
		aes(x = gen.change, ymin = count_lb, ymax = count_ub), 
		width = 1, 
		col = "black", linetype = "solid", size = 0.5, na.rm = TRUE
	)  +

	scale_fill_manual(
		name = "Virus sample",
		values = c(
			"B.1.36.16: 73NLt" = "#00A087FF", 
			"B.1.36.16: CV130" = "#91D1C2DD", 
			"AY.30: OTV54" = "#E64B35FF", 
			"AY.30: NH783" = "#F39B7FFF" 
		)
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
			"B.1.36.16: 73NLt" = NA, 
			"B.1.36.16: CV130" = NA, 
			"AY.30: OTV54" = NA, 
			"AY.30: NH783" = NA 
		)
	) +

	geom_text(
		data = data.frame(
			gen.change = "A.INS", count = 40,
			variant = factor("AY.30", levels(genetic_change_count_table$variant)),
			sample = factor("NH783", levels(genetic_change_count_table$sample)), 
			cell_line = factor("Calu-3", levels(genetic_change_count_table$cell_line)),
			exp_id = factor("AY.30: NH783")),
		label = "N/A",
		show.legend = FALSE
	) + 
	
	scale_x_discrete(
		name = "Genetic change", 
		limits = names(variant_table_dat %>% select(A.T:DEL.G) %>% select(-contains("INS.")) ),
		labels = setNames(
			gsub(
				"\\.",
				">",
				names(variant_table_dat %>% select(A.T:DEL.G) %>% select(-contains("INS.")) )
			),
			names(variant_table_dat %>% select(A.T:DEL.G) %>% select(-contains("INS.")))
		)
	) +
	scale_y_continuous(name = "Polymorphic site count") +
	facet_nested(
		cell_line ~ variant + sample, 
		labeller = labeller(cell_line = c("Clinical sample" = "Clinical\nsample", "Vero E6" = "Vero E6", "Vero E6-TMPRSS2" = "Vero E6/\nTMPRSS2", "Calu-3" = "Calu-3"))
	) +
	theme_light() + 
	theme(
		legend.position = "none",

		axis.text.x = element_text(angle = 90, vjust = 0, hjust=1, size = 6), axis.text.y = element_text(size = 6),

		strip.background =element_rect(fill = "white", color = "black"),
		strip.text = element_text(size = 8, color = "black", face = "bold")
	)

####################
#Make a composit plot 
####################
aligned_plots <- align_plots(
	var_count_through_passages_plot,
	genetic_change_count_plot,
	align = "v", axis = "lr"
)

Fig3_plot <- plot_grid(
	aligned_plots[[1]],
	aligned_plots[[2]],
	labels = c("a)", "b)"),
	ncol = 1, rel_heights = c(2, 4)
	)+
	theme(plot.background = element_rect(fill = "white", color = NA),  panel.border = element_blank())

ggsave(Fig3_genetic_change_count.png_filename,
	plot = Fig3_plot,
	width = 16, height = 16, units = "cm",
	dpi = 300)

ggsave(Fig3_genetic_change_count.svg_filename,
	plot = Fig3_plot,
	width = 16, height = 16, units = "cm",
	dpi = 300)

####################
#Comparing polymorphic site count distributions by likelihood ratio tests ..
####################
#subset the data to include only variant site counts from cultured samples
genetic_change_count_no_clinical_sample <- genetic_change_count_table %>% 
	filter(cell_line %in% c("Vero E6", "Vero E6-TMPRSS2", "Calu-3")) %>%
	select(variant:G.DEL, DEL.A:DEL.G) %>%
	gather(type, count, A.T:DEL.G)

#Construct models describing polymorphic site count distributions 
mdl_1 <- glm(formula = count ~ type + sample * cell_line , data = genetic_change_count_no_clinical_sample, family = "poisson") 
mdl_2 <- glm(formula = count ~ type + sample + cell_line , data = genetic_change_count_no_clinical_sample, family = "poisson") 
mdl_3 <- glm(formula = count ~ type, data = genetic_change_count_no_clinical_sample, family = "poisson") 

#Likelihood ratio tests
sink(Misc_info_var_site_count_by_cellline_filename, append = T)
cat("Fitting Poisson regression models to the count distributions\n")
cat("=========================\n")
cat(sprintf("model: %s\n", format(mdl_1$formula)))
summary(mdl_1)
cat("\n")

cat(sprintf("model: %s\n", format(mdl_2$formula)))
summary(mdl_2)
cat("\n")

cat(sprintf("model: %s\n", format(mdl_3$formula)))
summary(mdl_3)
cat("\n")

cat("Model comparisons\n")
cat("=========================\n")
cat("**Count distributions were significantly different among viral samples and cell culture conditions...**\n")
anova(mdl_3, mdl_2, test = "LRT")
cat("\n")
cat("**...and the effects were not simply additive**\n")
anova(mdl_2, mdl_1, test = "LRT")
cat("\n")
sink()

####################
#Comparing polymorphic site count distributions by likelihood ratio tests, considering base substitutions only
####################
#subset the data to include only variant site counts with base sub types found in cultured samples
genetic_change_count_no_clinical_sample_base_sub_only <- genetic_change_count_table %>% 
	filter(cell_line %in% c("Vero E6", "Vero E6-TMPRSS2", "Calu-3")) %>%
	select(variant:G.C) %>%
	gather(type, count, A.T:G.C)

#Construct the models describing polymorphic site count distributions 
mdl_1 <- glm(formula = count ~ type + sample * cell_line , data = genetic_change_count_no_clinical_sample_base_sub_only, family = "poisson") 
mdl_2 <- glm(formula = count ~ type + sample + cell_line , data = genetic_change_count_no_clinical_sample_base_sub_only, family = "poisson") 
mdl_3 <- glm(formula = count ~ type, data = genetic_change_count_no_clinical_sample_base_sub_only, family = "poisson") 

#Likelihood ratio tests
sink(Misc_info_var_site_count_by_cellline_filename, append = T)
cat("Fitting Poisson regression models to the count distributions, considering base substitutions only\n")
cat("=========================\n")
cat(sprintf("model: %s\n", format(mdl_1$formula)))
summary(mdl_1)
cat("\n")

cat(sprintf("model: %s\n", format(mdl_2$formula)))
summary(mdl_2)
cat("\n")

cat(sprintf("model: %s\n", format(mdl_3$formula)))
summary(mdl_3)
cat("\n")

cat("Model comparisons\n")
cat("=========================\n")
cat("**Count distributions were significantly different among viral samples and cell culture conditions...**\n")
anova(mdl_3, mdl_2, test = "LRT")
cat("\n")
cat("**...and the effects were not simply additive**\n")
anova(mdl_2, mdl_1, test = "LRT")
cat("\n")
sink()

####################
#Comparing polymorphic site freq distributions by likelihood ratio tests, considering base substitutions only
####################
genetic_change_freq_no_clinical_sample_base_sub_only <- 
genetic_change_count_no_clinical_sample_base_sub_only %>%
group_by(variant, sample, cell_line, passage, replicate) %>%
mutate(freq = count/sum(count)) %>% 
ungroup() %>% as.data.frame %>%
arrange(variant, sample, cell_line, passage, replicate) 

mdl_1 <- glm(formula = freq ~ type + sample * cell_line , data = genetic_change_freq_no_clinical_sample_base_sub_only, family= "quasibinomial") 
mdl_2 <- glm(formula = freq ~ type + sample + cell_line , data = genetic_change_freq_no_clinical_sample_base_sub_only, family = "quasibinomial") 
mdl_3 <- glm(formula = freq ~ type, data = genetic_change_freq_no_clinical_sample_base_sub_only, family = "quasibinomial") 

sink(Misc_info_var_site_count_by_cellline_filename, append = T)
cat("\n")
cat("Fitting quasibinomial regression models to the freq distributions, considering base substitutions only\n")
cat("=========================\n")
cat(sprintf("model: %s\n", format(mdl_1$formula)))
summary(mdl_1)
cat("\n")

cat(sprintf("model: %s\n", format(mdl_2$formula)))
summary(mdl_2)
cat("\n")

cat(sprintf("model: %s\n", format(mdl_3$formula)))
summary(mdl_3)
cat("\n")

cat("Model comparisons\n")
cat("=========================\n")
cat("**freq distributions were insignificantly different among viral samples and cell culture conditions...**\n")
anova(mdl_3, mdl_2, test = "LRT")
sink()


####################
#Set working directory
####################
path_to_wd = "path/to/code/dir" # path to the "code" folder
setwd(path_to_wd)

####################
#Load required libraries
####################
library(dplyr)
library(hash)
library(tidyverse)

library(ape)

library(ggplot2)
library(ggh4x) #facet_nested function
library(cowplot)

source("miscellaneous funcs.R")

####################
#Set paths to file
####################
#inputs
ct_value_dat_filename <- "data/ct values.txt"
virus_titre_dat_filename <- "data/infectious virus titres.txt"

#outputs
formated_ct_value_dat_filename <- "results/Supplementary Table 1.ct_value_by_passage_stock.txt"
Fig1_ct_and_virus_titre.raw.png_filename <- "results/Fig1_ct_and_virus_titre.raw.png"
Fig1_ct_and_virus_titre.raw.svg_filename <- "results/Fig1_ct_and_virus_titre.raw.svg"

####################
#Make Ct value plot
####################
#load ct value data
ct_value_dat <- read.table(ct_value_dat_filename, header = T, sep = "\t") %>% 
	set_factors #from "miscellaneous funcs.R"

#format the ct value table and write it to file
ct_value_dat %>% select(-variant_sample) %>% spread(passage, ct) %>%
	write.table(file = formated_ct_value_dat_filename, quote = F, sep = "\t", row.name = F)

#make the plot
ct_value_plot <- ct_value_dat %>% mutate(exp_id = paste(sample, cell_line, replicate, sep = "_")) %>%
	ggplot(aes(x = passage, y = ct, group = exp_id, col = variant_sample, linetype = replicate), shape = 16) +
	geom_point(size = 2) + 
	geom_line(size = 0.5) + 
	scale_color_manual(
		name = "Virus sample",
		values = c(
			"B.1.36.16: 73NLt" = "#00A087FF",
			"B.1.36.16: CV130" = "#91D1C2DD",
			"AY.30: OTV54" = "#E64B35FF",
			"AY.30: NH783" = "#F39B7FFF"
		),
		guide = guide_legend(order = 2)
	) +
	scale_linetype_manual(
		name = "Replicate",
		values = c( 
			"A" = "solid", 
			"B" = "dashed"
		), 
		guide = guide_legend(direction = "horizontal", order = 1)
	) +
	facet_grid(. ~ cell_line, 
		scales = "free", space = "free", 
		labeller = labeller(cell_line = c("Clinical sample" = "Clinical\nsample", "Vero E6" = "Vero E6", "Vero E6-TMPRSS2" = "Vero E6/TMPRSS2", "Calu-3" = "Calu-3"))) + 
	scale_x_discrete(name = "Passage") + 
	scale_y_continuous(name = "Ct value") +

	theme_light() + 
	theme(
		legend.position = "right",
		legend.spacing.y = unit(0, "lines"),
      	legend.title = element_text(size = 10),
	      legend.text = element_text(size = 8),

		strip.background = element_rect(fill = "white", color = "black"),
		strip.text = element_text(colour = "black", face = "bold"),

		panel.grid.minor = element_blank()
	)

####################
#Make infectious virus titre plot
####################
#load virus titre data
virus_titre_dat <- read.table(virus_titre_dat_filename, header = T, sep = "\t") %>% 
	as.data.frame() %>% 
	gather(passage, titre, P1:P4) %>%
	mutate(
		variant = factor(variant, levels = c("B.1.36.16", "AY.30")),
		sample = factor(sample, levels = c("73NLt", "CV130", "OTV54", "NH783")),
		cell_line = factor(cell_line, levels = c("Vero E6", "Vero E6-TMPRSS2", "Calu-3")),
		passage = factor(passage, levels = c("P1", "P2", "P3", "P4"), ordered = TRUE),
		t = as.numeric(passage),
		variant_sample = factor(paste(variant, sample, sep = ": "), levels = c("B.1.36.16: 73NLt", "B.1.36.16: CV130", "AY.30: OTV54", "AY.30: NH783")),
		titre = as.numeric(titre)
	)

#linear model fittings, log(P2-P4 titer) ~ t
mdl_pred <- NULL
mdl_coef <- NULL
for (cell_line_id in levels(virus_titre_dat$cell_line)){
	for (sample_id in levels(virus_titre_dat$sample)){
		if (cell_line_id == "Calu-3" & (sample_id == "OTV54" | sample_id == "NH783")){
			print("do not thing")
		} else {
			virus_titre_dat_pred <- virus_titre_dat %>%
				filter(t>1) %>%
				filter(sample == sample_id, cell_line == cell_line_id)
			
			mdl <- lm(formula = log(titre) ~ t, data = virus_titre_dat_pred) 

			mdl_pred <- rbind(
				mdl_pred, 
				cbind(
					virus_titre_dat_pred, 
					pred = exp(predict(mdl, virus_titre_dat_pred, interval = "confidence"))
				)
			)
			mdl_coef <- rbind(
				mdl_coef, 
				cbind(
					virus_titre_dat_pred %>% select(-t) %>% unique(),
					intercept = mdl$coef[1],
					t = mdl$coef[2],
					r2 = summary(mdl)$r.squared,
					p = summary(mdl)$coef[2,4]
				)
			)
		}
	}
}
mdl_pred <- mdl_pred %>% 
select(variant, sample, cell_line, t, variant_sample, pred.fit, pred.lwr, pred.upr) %>%
distinct

#make the plot
eq_text <- mdl_coef %>% 
	mutate(
		intercept = round(intercept, 2), 
		slope = ifelse(
			t > 0, 
			sprintf(" + %s", round(t, 2)),
			sprintf(" - %s", round(abs(t), 2))
		),
		r2 = round(r2, 2),

		sig = ifelse(
			p < 0.0001,
			"***",
			ifelse(
				p < 0.001,
				"**",
				ifelse(
					p < 0.005,
					"*",
					""
				)
			)
		),
		p = ifelse(
			p < 0.0001,
			"p < 0.0001",
			ifelse(
				p < 0.001,
				"p < 0.001",
				ifelse(
					p < 0.005,
					"p < 0.005",
					sprintf("p = %s", round(p, 4))
				)
			)
		)
	) %>%
	select(variant, sample, cell_line, variant_sample, intercept, slope, r2, p, sig) %>% unique() %>%
	mutate(label = paste("ln(titre) = ", intercept, slope, " × passage", ";\n", expression("r^2 = "), r2, "; ", p, sig, sep = "") )

inf_virus_titre_plot <- 
	ggplot(
		data = virus_titre_dat,
		aes(x = t, y = titre, col = variant_sample, shape = variant_sample)
	) +
	geom_point(size = 2) + 

	scale_color_manual(
		name = "Virus sample",
		values = c(
			"B.1.36.16: 73NLt" = "#00A087FF",
			"B.1.36.16: CV130" = "#91D1C2DD",
			"AY.30: OTV54" = "#E64B35FF",
			"AY.30: NH783" = "#F39B7FFF"
		)
	) +
	scale_shape_manual(
		name = "Virus sample",
		values = c(
			"B.1.36.16: 73NLt" = 16, 
			"B.1.36.16: CV130" = 16,
			"AY.30: OTV54" = 16, 
			"AY.30: NH783" = 16
		)
	) +
	geom_line(
		data = mdl_pred,
		aes(y = pred.fit, linetype = variant_sample),
		size = 0.5
	) +

	geom_errorbar(
		data = mdl_pred,
		aes(	y = pred.fit,
			ymin = pred.lwr,
			ymax = pred.upr,
			col = variant_sample
		), 
		width = .2
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

	geom_text(
		data = eq_text ,
		mapping = aes(x = 1, y = 10**2, label = label),
		hjust = 0, vjust = 1, 
		size = 2.5, color = "black"
		) +

	geom_text(
		data = data.frame(
			t = 2.5,
			titre = 10**3,
			variant = factor("AY.30", levels(virus_titre_dat$variant)),
			sample = factor("NH783", levels(virus_titre_dat$sample)),
			variant_sample = factor("AY.30: NH783", levels(virus_titre_dat$variant_sample)),
			cell_line = factor("Calu-3", levels(virus_titre_dat$cell_line))
		),
		label = "N/A",
		show.legend = FALSE
	) + 

	scale_x_continuous(name = "Passage", labels = c("1" = "P1", "2" = "P2", "3" = "P3", "4" = "P4")) +
	scale_y_continuous(name = "Infectious virus titre (FFU/mL)", trans = 'log10', limits = c(1,9E5), expand = c(0, 0)) + 
	facet_nested(
		variant + sample ~ cell_line,
		labeller = labeller(cell_line = setNames(c("Vero E6", "Vero E6/TMPRSS2", "Calu-3"), levels(virus_titre_dat$cell_line) ) )
	) + 

	theme_light() + 
	theme(
		legend.position = "none",
		strip.background = element_rect(fill = "white", color = "black"),
		strip.text = element_text(colour = "black", face = "bold"),
	)

####################
#Make a composit plot of Ct value, and infectious virus titre
####################
aligned_left_plots <- align_plots(
	ggplot() + theme_void(), # a blank plot for Fig 1A, Overview of SARS-CoV-2 propagation experiments
	ct_value_plot,
	inf_virus_titre_plot,
	align = "v", axis = "l"
)

Fig1_ct_and_virus_titre <- plot_grid(
	aligned_left_plots[[1]],
	aligned_left_plots[[2]],
	aligned_left_plots[[3]],
	labels = c("a)", "b)", "c)"),
	ncol = 1, rel_heights = c(1, 1, 2.5)
) +
theme(plot.background = element_rect(fill = "white", color = NA),  panel.border = element_blank())

ggsave(Fig1_ct_and_virus_titre.raw.png_filename, # Fig 1A made elsewhere; added to the plot manually
	plot = Fig1_ct_and_virus_titre,
	width = 16*1.25, height = 15*1.25, units = "cm",
	dpi = 300)
ggsave(Fig1_ct_and_virus_titre.raw.svg_filename, # Fig 1A made elsewhere; added to the plot manually
	plot = Fig1_ct_and_virus_titre,
	width = 16*1.25, height = 15*1.25, units = "cm",
	dpi = 300)

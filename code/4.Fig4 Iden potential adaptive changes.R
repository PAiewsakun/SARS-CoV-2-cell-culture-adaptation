####################
#Set working directory
####################
path_to_wd = "path/to/code/dir" # path to the "code" folder
setwd(path_to_wd)

####################
#Load required libraries
####################
library(dplyr)
library(tibble) #rownames_to_column
library(forcats) #fct_recode
library(tidyr) #unite
library(R.utils) #insert function
library(stringr) #str_replace_all function

library(ape)
library(lme4qtl) #relmatGlmer

library(ggplot2)
library(ggh4x) #facet_nested function
library(ggpmisc) #geom_table
library(ggnewscale) #new_scale
library(cowplot) #align_plots

source("miscellaneous funcs.R")

####################
#Set paths to file
####################
#inputs
variant_table_dat_filename <- "data/variant table.txt"
ref_annotations_filename <- "data/ref genome annotation.txt"
clinical_sample_ML_tree_filename <- "data/clinical viral sample ML tree.treefile"

#outputs
cand_adap_site_info_filename <- "results/Supplementary Data 3.cand_adap_site_info.raw.required formatting.txt"
Fig4_mut_adv_coef_by_pos_plot.png_filename <- "results/Fig4_mut_adv_coef_by_pos_plot.png"
Fig4_mut_adv_coef_by_pos_plot.svg_filename <- "results/Fig4_mut_adv_coef_by_pos_plot.svg"

Ext_Dat_FigS3_temp_dyn_of_adap_changes.E6.png_filename <- "results/Ext Dat FigS3_temp_dyn_of_adap_changes.E6.png"
Ext_Dat_FigS3_temp_dyn_of_adap_changes.E6.svg_filename <- "results/Ext Dat FigS3_temp_dyn_of_adap_changes.E6.svg"
Ext_Dat_FigS4_temp_dyn_of_adap_changes.TM.png_filename <- "results/Ext Dat FigS4_temp_dyn_of_adap_changes.TM.png"
Ext_Dat_FigS4_temp_dyn_of_adap_changes.TM.svg_filename <- "results/Ext Dat FigS4_temp_dyn_of_adap_changes.TM.svg"
Ext_Dat_FigS5_temp_dyn_of_adap_changes.C3.png_filename <- "results/Ext Dat FigS5_temp_dyn_of_adap_changes.C3.png"
Ext_Dat_FigS5_temp_dyn_of_adap_changes.C3.svg_filename <- "results/Ext Dat FigS5_temp_dyn_of_adap_changes.C3.svg"

####################
#Set parameter values
####################
#set depht, mutation frequency, mut_adv_coef, and p-value thresholds
depth_lw_threshold <- 30
fMut_lw_threshold <- 0.05
mut_adv_coef_lw_threshold <- 0
adj_p_threshold <- 0.05/1/29644/3 #Bonferroni multiple-testing corrected p value threshold from 0.05 to 0.05 / 29,644 sites / 3 cell lines = 5.62×10-7%

####################
#Compute mutant variant frequencies
####################
#load variant data
variant_table_dat <- read.table(variant_table_dat_filename, header = T, sep = "\t") %>% 
	set_factors %>% #from "miscellaneous funcs.R"
	mutate(
		t = as.numeric(passage) - 1,

		fA = A/depth,
		fT = T/depth,
		fC = C/depth,
		fG = G/depth,
		fINS = INS/depth,
		fDEL = DEL/depth
	) %>% replace(is.na(.), 0) %>% 
	arrange(position, variant, sample, cell_line, replicate, passage) 

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

#compute mutant variant frequencies
variant_table_dat$fMut <- 1 - variant_table_dat$fOri.var #"The major variant present in the original clinical sample was taken as the original variant, and any other variants that were not the original variant were collectively grouped together as the mutant variant..."

####################
#Perform mixed effects logistic regressions 
####################
#for characterising mutations
ref_seq <- toupper(as.character(read.GenBank("NC_045512.2"))[[1]])
ref_annotations <- read.table(ref_annotations_filename, header = TRUE, sep = "\t", quote = "")

variant_table_dat_for_plotting <- NULL
cand_adap_site_info <- NULL
model_fitting_err_sites <- NULL

for (cell_line_id in c("Vero E6", "Vero E6-TMPRSS2", "Calu-3")){
	#Compute variance-covariance (vcv) matrix among clinical viral samples
	if (cell_line_id == "Calu-3"){
		#Load the ML tree
		phy <- read.tree(file = clinical_sample_ML_tree_filename)
		phy <- root(phy, outgroup = "NC_045512.2", resolve.root = TRUE)
		phy <- drop.tip(phy, c("NC_045512.2", "NH783")) #NH783 could not be propagated on Calu-3 cells..

		#Compute phy vcv matrix
		Vphy <- vcv(phy)
		Vphy <- Vphy/det(Vphy)^(1/nrow(Vphy)) #normalise the vcv matrix to have a determinant of 1
	} else {
		#Load the ML tree
		phy <- read.tree(file = clinical_sample_ML_tree_filename)
		phy <- root(phy, outgroup = "NC_045512.2", resolve.root = TRUE)
		phy <- drop.tip(phy, c("NC_045512.2"))

		#Compute phy vcv matrix
		Vphy <- vcv(phy)
		Vphy <- Vphy/det(Vphy)^(1/nrow(Vphy)) #normalise the vcv matrix to have a determinant of 1
	}
	
	#get sites with fMut >= fMut_lw_threshold in at least one dataset,...
	#...given that corresponding sequencing depth >= depth_lw_threshold
	position_to_be_analysed <- variant_table_dat %>% filter(cell_line %in% c(cell_line_id, "Clinical sample")) %>% 
		mutate(fMut = ifelse(depth >= depth_lw_threshold, fMut, 0)) %>%
		group_by(position) %>%
		summarise(
			position_to_be_analysed = any(fMut >= fMut_lw_threshold)
		) %>%
		ungroup() %>% as.data.frame() %>% 
		filter(position_to_be_analysed == TRUE) %>% pull(position)

	for (position_id in position_to_be_analysed){
		print(sprintf("Cell line: %s; position: %s", cell_line_id, position_id))

		#subset 'variant_table_dat' to include only those from 'cell_line_id' at 'position_id'
		subset_variant_table_dat <- variant_table_dat %>% filter(cell_line == cell_line_id, position == position_id)

		#duplicate original clinical data for each experimental replicate
		subset_variant_table_dat <- rbind(
			subset_variant_table_dat,
			left_join(
				subset_variant_table_dat %>% select(sample, replicate) %>% unique, 
				variant_table_dat %>% filter(cell_line == "Clinical sample", position == position_id) %>% select(-replicate) %>% mutate(cell_line = cell_line_id), 
				by = c("sample" = "sample")
			) %>%
			relocate(sample, .after = seq_name) %>% relocate(replicate, .after = passage) 
		) %>% arrange(variant, sample, replicate, passage)

		#set mu_freq with 0 depth to NA
		subset_variant_table_dat[subset_variant_table_dat$depth < depth_lw_threshold, "fMut"] <- NA

		#make mu_nt and mu_aa coloumn for recording mutations
		subset_variant_table_dat <- subset_variant_table_dat %>% mutate(mu_nt_aa = "")
		
		#mixed model fitting with lmer4qlt
		err <- 0
		M0 <- NULL #the initial mutation frequencies (i.e. the model intercepts) to vary among the B.1.36.16 and AY.30 variants, while assuming that the mutant variant bares no selective advantage over the original variant.
		M1 <- NULL #the initial mutation frequencies and the ?s values were estimated from the data, and were allowed to vary among the B.1.36.16 and AY.30 variants. The ?s values were assumed to be shared and do not randomly vary significantly among viral samples of the same variant and experimental replicates.
		M2 <- NULL #the initial mutation frequencies and the ?s values were estimated from the data, and were allowed to vary among the B.1.36.16 and AY.30 variants. The ?s values were allowed to vary randomly among viral samples and experimental replicates.

		M0VsM1_anova_res <- NULL
		M1VsM2_anova_res <- NULL

		best_fit_mdl_name <- NULL
		best_fit_mdl <- NULL

	tryCatch({ ######### try model fitting. There will be some error since some sites have very low sequencing depths #########
		M0 <- relmatGlmer(
			fMut ~ variant - 1 + (1|sample/replicate), 
			data = subset_variant_table_dat, 
			weights = depth,
			relmat = list(sample = Vphy),
			family = binomial
		)

		M1 <- relmatGlmer(
			fMut ~ variant - 1 + variant:t + (1|sample/replicate), 
			data = subset_variant_table_dat, 
			weights = depth,
			relmat = list(sample = Vphy),
			family = binomial
		)

		M2 <- relmatGlmer(
			fMut ~ variant - 1 + variant:t + (1+t|sample/replicate), 
			data = subset_variant_table_dat, 
			weights = depth,
			relmat = list(sample = Vphy),
			family = binomial
		)

		# Model comparisons by ANOVA test
		M0VsM1_anova_res <- anova(M0, M1)
		M0VsM1_anova_chisq <- M0VsM1_anova_res$"Chisq"[2]
		M0VsM1_anova_df <- M0VsM1_anova_res$"Df"[2]
		M0VsM1_anova_p <- M0VsM1_anova_res$"Pr(>Chisq)"[2]

		M1VsM2_anova_res <- anova(M1, M2)
		M1VsM2_anova_chisq <- M1VsM2_anova_res$"Chisq"[2]
		M1VsM2_anova_df <- M1VsM2_anova_res$"Df"[2]
		M1VsM2_anova_p <- M1VsM2_anova_res$"Pr(>Chisq)"[2]

		best_fit_mdl_name <- ifelse(M0VsM1_anova_p > adj_p_threshold, "M0", ifelse(M1VsM2_anova_p > adj_p_threshold, "M1", "M2"))
		if (best_fit_mdl_name == "M1") {best_fit_mdl = M1}
		if (best_fit_mdl_name == "M2") {best_fit_mdl = M2}

		#make plot and keep the results when best-fit model == M1|M2 and ...
		#...selection coeffecients of at least one virus variant are > 0 and ...
		#...the correposnding variant with a positive selection coeffecient has fMut >= fMut_lw_threshold in at least one dataset
		if(
			best_fit_mdl_name %in% c("M1", "M2") & 
			(
				coef(summary(best_fit_mdl))["variantB.1.36.16:t", "Estimate"] > mut_adv_coef_lw_threshold & any(subset_variant_table_dat %>% filter(variant == "B.1.36.16") %>% pull(fMut) >= fMut_lw_threshold, na.rm = T) |
				coef(summary(best_fit_mdl))["variantAY.30:t", "Estimate"] > mut_adv_coef_lw_threshold & any(subset_variant_table_dat %>% filter(variant == "AY.30") %>% pull(fMut) >= fMut_lw_threshold, na.rm = T)
			) 
		){
			#print best fit model
			print(summary(best_fit_mdl))

			#get mut_adv_coef and associated info for the two variants
			mut_adv_coef_B.1.36.16 <- coef(summary(best_fit_mdl))["variantB.1.36.16:t", "Estimate"]
			mut_adv_coef_se_B.1.36.16 <- coef(summary(best_fit_mdl))["variantB.1.36.16:t", "Std. Error"]
			mut_adv_coef_z_B.1.36.16 <- coef(summary(best_fit_mdl))["variantB.1.36.16:t", "z value"]
			mut_adv_coef_p_B.1.36.16 <- coef(summary(best_fit_mdl))["variantB.1.36.16:t", "Pr(>|z|)"]

			mut_adv_coef_AY.30 <- coef(summary(best_fit_mdl))["variantAY.30:t", "Estimate"]
			mut_adv_coef_se_AY.30 <- coef(summary(best_fit_mdl))["variantAY.30:t", "Std. Error"]
			mut_adv_coef_z_AY.30 <- coef(summary(best_fit_mdl))["variantAY.30:t", "z value"]
			mut_adv_coef_p_AY.30 <- coef(summary(best_fit_mdl))["variantAY.30:t", "Pr(>|z|)"]

			#get nucleotide and amino acid changes
			for(j in 1:nrow(subset_variant_table_dat)){
				ori.var <- subset_variant_table_dat[j, "ori.var"]
				fMut_profile <- subset_variant_table_dat[j, c("fA", "fT", "fC", "fG", "fINS", "fDEL")]
				fMut_profile <- fMut_profile[! names(fMut_profile) %in% c(paste("f", ori.var, sep = ""))]
				mu.vars <- gsub("f", "", names(fMut_profile)[fMut_profile > fMut_lw_threshold])

				#substitution type
				sub.mu.vars <- mu.vars[!mu.vars %in% c("DEL", "INS")]
				if(length(sub.mu.vars) != 0){
					subset_variant_table_dat[j, "mu_nt_aa"] <- 
						paste(
							paste(
								paste(prettyNum(position_id, big.mark = ","), ori.var, ">", sub.mu.vars, sep = ""),
								"(", 
								nucleotide_mu_to_protein_mu(
									position = position_id, 
									ori = ori.var, 
									mu = sub.mu.vars,
									ref = ref_seq,
									ref_annotations = ref_annotations %>% filter(type == "gene")),
								")",
								sep = ""
							), collapse = "; "
						)
				}

				#del type
				if("DEL" %in% mu.vars){
					subset_variant_table_dat[j, "mu_nt_aa"] <- 
						paste(
							remove_empty_strings_from_a_vector(
								c(
									subset_variant_table_dat[j, "mu_nt_aa"],
									paste(prettyNum(position_id, big.mark = ","), "del", ori.var, sep = "")
								)
							),
							collapse = "; "
						)
				}

				#ins type 
				if("INS" %in% mu.vars){
					subset_variant_table_dat[j, "mu_nt_aa"] <- 
						paste(
							remove_empty_strings_from_a_vector(
								c(
									subset_variant_table_dat[j, "mu_nt_aa"],
									paste(prettyNum(position_id, big.mark = ","), "_", prettyNum(position_id + 1, big.mark = ","), "ins", "XXX", sep = "") #require manual curation!!!!
								)
							),
							collapse = "; "
						)
				}

				#set mu_nt and mu_aa to "", if..
				if(
					subset_variant_table_dat[j, "depth"] < depth_lw_threshold |
					(subset_variant_table_dat[j, "variant"] == "B.1.36.16" & mut_adv_coef_B.1.36.16 < mut_adv_coef_lw_threshold) |
					(subset_variant_table_dat[j, "variant"] == "AY.30" & mut_adv_coef_AY.30 < mut_adv_coef_lw_threshold)
				){
					subset_variant_table_dat[j, "mu_nt_aa"] <- ""
				}
			}

			variant_table_dat_for_plotting <- rbind(
				variant_table_dat_for_plotting,
				subset_variant_table_dat %>% 
					select(variant, sample, cell_line, passage, replicate, position, fMut, mu_nt_aa) %>%
					arrange(variant, sample, cell_line, replicate, passage)
			)

			#keep information regarding candidate adaptive sites
			mu_nt_aa_B.1.36.16 <- paste(subset_variant_table_dat %>% filter(variant == "B.1.36.16") %>% pull(mu_nt_aa) %>% remove_empty_strings_from_a_vector %>% sapply(function(x) strsplit(x, "; ")) %>% unlist %>% unique, collapse = "; ")
			mu_nt_aa_AY.30 <- paste(subset_variant_table_dat %>% filter(variant == "AY.30") %>% pull(mu_nt_aa) %>% remove_empty_strings_from_a_vector %>% sapply(function(x) strsplit(x, "; ")) %>% unlist %>% unique, collapse = "; ")

			gene <- get_gene_from_pos(position_id, ref_annotations %>% filter(type == "gene"))
			gene <- ifelse(is.na(gene), "", gene)
			subdomain <- get_region_from_pos(position_id, ref_annotations %>% filter(type == "subregion"))
			subdomain <- ifelse(is.na(subdomain), "", subdomain)
			
			cand_adap_site_info <- rbind(
				cand_adap_site_info, 
				c(
					c(
						"cell_line" = cell_line_id, 
						"gene" = gene,
						"subdomain" = subdomain,
						"position" = position_id, 

						"M0VsM1_anova_chisq" = M0VsM1_anova_chisq,
						"M0VsM1_anova_df" = M0VsM1_anova_df,
						"M0VsM1_anova_p" = M0VsM1_anova_p,

						"M1VsM2_anova_chisq" = M1VsM2_anova_chisq,
						"M1VsM2_anova_df" = M1VsM2_anova_df,
						"M1VsM2_anova_p" = M1VsM2_anova_p,

						"best_fit_mdl_name" = best_fit_mdl_name,

						"mut_adv_coef_B.1.36.16" = mut_adv_coef_B.1.36.16,
						"mut_adv_coef_se_B.1.36.16" = mut_adv_coef_se_B.1.36.16,
						"mut_adv_coef_z_B.1.36.16" = mut_adv_coef_z_B.1.36.16,
						"mut_adv_coef_p_B.1.36.16" = mut_adv_coef_p_B.1.36.16,

						"mut_adv_coef_AY.30" = mut_adv_coef_AY.30,
						"mut_adv_coef_se_AY.30" = mut_adv_coef_se_AY.30,
						"mut_adv_coef_z_AY.30" = mut_adv_coef_z_AY.30,
						"mut_adv_coef_p_AY.30" = mut_adv_coef_p_AY.30,

						"adap_site_cand_B.1.36.16" = mut_adv_coef_B.1.36.16 > mut_adv_coef_lw_threshold & any(subset_variant_table_dat %>% filter(variant == "B.1.36.16") %>% pull(fMut) >= fMut_lw_threshold, na.rm = T),
						"adap_site_cand_AY.30" = mut_adv_coef_AY.30 > mut_adv_coef_lw_threshold & any(subset_variant_table_dat %>% filter(variant == "AY.30") %>% pull(fMut) >= fMut_lw_threshold, na.rm = T),

						"mu_nt_aa_B.1.36.16" = mu_nt_aa_B.1.36.16,
						"mu_nt_aa_AY.30" = mu_nt_aa_AY.30
					)
				)
			)
		}
	}, error = function(e){
		err <<- 1
	})
	if(err == 1){
		model_fitting_err_sites[[cell_line_id]] <- c(model_fitting_err_sites[[cell_line_id]], position_id)
		err <- 0
		next
	}
	}
}

#fotmat the df cand_adap_site_info a bit
cand_adap_site_info <- cand_adap_site_info %>% as.data.frame() %>%
	mutate(
		cell_line = factor(cell_line, levels = c("Vero E6", "Vero E6-TMPRSS2", "Calu-3")),
		position = as.numeric(position),

		M0VsM1_anova_chisq = as.numeric(M0VsM1_anova_chisq),
		M0VsM1_anova_df = as.integer(M0VsM1_anova_df),
		M0VsM1_anova_p = as.numeric(M0VsM1_anova_p),

		M1VsM2_anova_chisq = as.numeric(M1VsM2_anova_chisq),
		M1VsM2_anova_df = as.integer(M1VsM2_anova_df),
		M1VsM2_anova_p = as.numeric(M1VsM2_anova_p),

		best_fit_mdl_name = factor(best_fit_mdl_name, levels = c("M1", "M2")),

		mut_adv_coef_B.1.36.16 = as.numeric(mut_adv_coef_B.1.36.16),
		mut_adv_coef_se_B.1.36.16 = as.numeric(mut_adv_coef_se_B.1.36.16),
		mut_adv_coef_z_B.1.36.16 = as.numeric(mut_adv_coef_z_B.1.36.16),
		mut_adv_coef_p_B.1.36.16 = as.numeric(mut_adv_coef_p_B.1.36.16),

		mut_adv_coef_AY.30 = as.numeric(mut_adv_coef_AY.30),
		mut_adv_coef_se_AY.30 = as.numeric(mut_adv_coef_se_AY.30),
		mut_adv_coef_z_AY.30 = as.numeric(mut_adv_coef_z_AY.30),
		mut_adv_coef_p_AY.30 = as.numeric(mut_adv_coef_p_AY.30),

		adap_site_cand_B.1.36.16 = as.logical(adap_site_cand_B.1.36.16),
		adap_site_cand_AY.30 = as.logical(adap_site_cand_AY.30)
	)

#data pre-processing
#add alpha and positive_exp_id to the df
variant_table_dat_for_plotting <- variant_table_dat_for_plotting %>% 
	mutate(exp_id_pos = paste(sample, cell_line, replicate, position, sep = "_"))
variant_table_dat_for_plotting <- left_join(
	variant_table_dat_for_plotting,
	variant_table_dat_for_plotting %>% 
		group_by(exp_id_pos) %>% 
		summarise(
			alpha = ifelse(all(mu_nt_aa == ""), 0.50, 1),
			positive_exp_id = ifelse(all(mu_nt_aa == ""), "", paste(sample, replicate, sep = ":")) 
		) %>% 
		ungroup() %>% as.data.frame(),
	by = c("exp_id_pos"))

#dirty, quick, manual fix for position 22419!!!!!!; total fMut >= fMut_lw_threshold, but no individual variant has freq >= fMut_lw_threshold 
variant_table_dat_for_plotting[
	variant_table_dat_for_plotting$exp_id_pos == "CV130_Calu-3_B_24419" & variant_table_dat_for_plotting$passage == "P4", 
	c("mu_nt_aa", "positive_exp_id")
	] <- c("24,419A>C(N953H); 24,419A>G(N953D); 24,419A>T(N953Y)", "CV130:B")

variant_table_dat_for_plotting[variant_table_dat_for_plotting$exp_id_pos == "CV130_Calu-3_B_24419", "alpha"] <- 1.00

cand_adap_site_info[
	cand_adap_site_info$cell_line == "Calu-3" & cand_adap_site_info$position == 24419,
	"mu_nt_aa_B.1.36.16"
	] <- "24,419A>C(N953H); 24,419A>G(N953D); 24,419A>T(N953Y)"

#dirty, quick, manual fix for position 22205!!!!!!; insertion type
variant_table_dat_for_plotting <- variant_table_dat_for_plotting %>% 
	mutate(mu_nt_aa = ifelse(mu_nt_aa == "22,205_22,206insXXX", "22,205_22,206insGA AGA GAT CGC CAT T(D215G_L216insRDRHY)", mu_nt_aa))

cand_adap_site_info <- cand_adap_site_info %>%
	mutate(mu_nt_aa_B.1.36.16 = ifelse(mu_nt_aa_B.1.36.16 == "22,205_22,206insXXX", "22,205_22,206insGA\nAGA GAT CGC CAT T\n(D215G_L216insRDRHY)", mu_nt_aa_B.1.36.16))

#add positive_exp_id_B.1.36.16 and positive_exp_id_AY.30 to cand_adap_site_info
cand_adap_site_info <- left_join(
	cand_adap_site_info,
	variant_table_dat_for_plotting %>% filter(variant == "B.1.36.16") %>%
	group_by(cell_line, position) %>%
	summarise(
		positive_exp_id_B.1.36.16 = {
			x <- unique(positive_exp_id[positive_exp_id != ""])
			n <- length(x)
			if(n == 0) {
				""
			}else if(n == 1) {
				x
			}else {
				x <- insert(x, seq(2,length(x),3), "\n")
				x <- insert(x, seq(1,length(x),1)[-1], ", ")
				x <- paste(x, collapse = "")
				gsub("\n, ", "\n", x)
			}
		}
	) %>%
	ungroup() %>% as.data.frame,
	by = c("cell_line", "position"))

cand_adap_site_info <- left_join(
	cand_adap_site_info,
	variant_table_dat_for_plotting %>% filter(variant == "AY.30") %>%
	group_by(cell_line, position) %>%
	summarise(
		positive_exp_id_AY.30 = {
			x <- unique(positive_exp_id[positive_exp_id != ""])
			n <- length(x)
			if(n == 0) {
				""
			}else if(n == 1) {
				x
			}else {
				x <- insert(x, seq(2,length(x),3), "\n")
				x <- insert(x, seq(1,length(x),1)[-1], ", ")
				x <- paste(x, collapse = "")
				gsub("\n, ", "\n", x)
			}
		}
	) %>%
	ungroup() %>% as.data.frame,
	by = c("cell_line", "position"))

#format mu_nt_aa for plotting
cand_adap_site_info <- cand_adap_site_info %>%
mutate(
	mu_nt_aa_B.1.36.16 = str_replace_all(mu_nt_aa_B.1.36.16, "; ", ";\n"),
	mu_nt_aa_AY.30 = str_replace_all(mu_nt_aa_AY.30, "; ", ";\n")
)

#write "cand_adap_site_info" to file
cand_adap_site_info %>% 
mutate(
	positive_exp_id_B.1.36.16 = gsub("\n", "", positive_exp_id_B.1.36.16),
	positive_exp_id_AY.30 = gsub("\n", "", positive_exp_id_AY.30),

	mu_nt_aa_B.1.36.16 = gsub("\n", " ", mu_nt_aa_B.1.36.16),
	mu_nt_aa_AY.30 = gsub("\n", " ", mu_nt_aa_AY.30),
) %>% 
write.table( 
	cand_adap_site_info_filename,
	row.name = F, quote = F, sep = "\t"
)

####################
#Ploting sites with potential adaptation mutations
####################
#make mut_adv_coef_by_pos_plot
mut_adv_coef_dat <- cand_adap_site_info %>%
	pivot_longer(
		cols = c(mut_adv_coef_B.1.36.16, mut_adv_coef_AY.30),
		names_to = "variant", 
		names_prefix = "mut_adv_coef_",
		values_to = "mut_adv_coef"
	) %>% 
	filter((variant == "B.1.36.16" & adap_site_cand_B.1.36.16)|(variant == "AY.30" & adap_site_cand_AY.30)) %>%
	mutate(
		cell_line = fct_recode(cell_line, "Vero E6/TMPRSS2" = "Vero E6-TMPRSS2"), 
		p = ifelse(best_fit_mdl_name == "M1", M0VsM1_anova_p, M1VsM2_anova_p),
		p = ifelse(p < 1e-30, 1e-30, p),
		mut_adv_coef = ifelse(mut_adv_coef > 2, 2, mut_adv_coef)
	) %>% 
	select(variant, cell_line, position, best_fit_mdl_name, p, mut_adv_coef) %>% as.data.frame 

mut_adv_coef_by_pos_plot <- mut_adv_coef_dat %>% 
	ggplot(aes(x = position, y = mut_adv_coef, fill = best_fit_mdl_name)) +
	geom_col(width = 20, alpha = 0.5) +

	scale_fill_manual(
		name = "\u0394s consistency",
		values = c(
			"M1" = "#3C5488FF",
			"M2" = "#DC0000FF"
		),
		labels = c("M1" = "Consistent \u0394s", "M2" = "Varying \u0394s")
	) +

	geom_hline(yintercept = 0) +

	geom_text(
		data = mut_adv_coef_dat %>% 
			group_by(cell_line, variant) %>%
			summarise(
				M1 = sum(best_fit_mdl_name == "M1"),
			) %>% ungroup() %>% as.data.frame %>%
			mutate(
				label = paste(
					"Consistent \u0394s: ", M1, " sites",
					sep = ""
				)
			),
		mapping = aes(x = 100, y = 2, label = label),
		col = "#3C5488FF",
		size = 1.75,
		hjust = 0, vjust = 1, 
	) + 

	geom_text(
		data = mut_adv_coef_dat %>% 
			group_by(cell_line, variant) %>%
			summarise(
				M2 = sum(best_fit_mdl_name == "M2"),
			) %>% ungroup() %>% as.data.frame %>%
			mutate(
				label = paste(
					"varying \u0394s: ", M2, " sites",
					sep = ""
				)
			),
		mapping = aes(x = 100, y = 1.5, label = label),
		col = "#DC0000FF",
		size = 1.75,
		hjust = 0, vjust = 1, 
	) + 

	facet_nested(cell_line + variant ~ .) +

	scale_x_continuous(name = "Position", limits = c(1,29903), expand = c(0, 0)) +
	scale_y_continuous(name = "Mutation selective advantage coeffecient (\u0394s)", limits = c(0, 2)) +
	guides(fill = guide_legend(title.position = "top")) +

	theme_classic() +
	theme(
		axis.line.x = element_blank(),

		axis.title = element_text(size = 8),
		axis.text = element_text(size = 6),
		strip.text = element_text(size = 6),

		legend.position = "bottom",
		legend.background = element_blank(),
		legend.margin = margin(t = 0, b = 0, l = 0, r = 0, unit = "cm"),
		legend.key.height = unit(0.25, 'cm'),
		legend.key.width = unit(0.25, 'cm'),
		legend.title = element_text(size = 6), 
		legend.text = element_text(size = 5)
	)

leg <- get_legend(mut_adv_coef_by_pos_plot)
mut_adv_coef_by_pos_plot <- mut_adv_coef_by_pos_plot + theme(legend.position = "none")

#Make a composit plot of genome structure, and sequencing depth
mut_adv_coef_by_pos_plot <- align_plots(
	make_ref_genome_structure_plot(text_size = 1.75),
	mut_adv_coef_by_pos_plot,
	align = "v", axis = "l"
)

mut_adv_coef_by_pos_plot <- align_plots(
	mut_adv_coef_by_pos_plot[[1]],
	mut_adv_coef_by_pos_plot[[2]],
	align = "v", axis = "lr"
)

mut_adv_coef_by_pos_plot <- ggdraw() + 
	draw_plot(mut_adv_coef_by_pos_plot[[2]],	x = 0.00, y = 0.00, width = 1.0, height = 0.90) +
	draw_plot(leg,					x = 0.80, y = 0.02, width = 0.1, height = 0.05) +
	draw_plot(mut_adv_coef_by_pos_plot[[1]],	x = 0.00, y = 0.875, width = 1.0, height = 0.125)

mut_adv_coef_by_pos_plot <- mut_adv_coef_by_pos_plot +
	theme(plot.background = element_rect(fill = "white", color = NA),  panel.border = element_blank())

ggsave(Fig4_mut_adv_coef_by_pos_plot.png_filename,
	plot = mut_adv_coef_by_pos_plot,
	width = 16, height = 9, units = "cm",
	dpi = 300)

ggsave(Fig4_mut_adv_coef_by_pos_plot.svg_filename,
	plot = mut_adv_coef_by_pos_plot,
	width = 16, height = 9, units = "cm",
	dpi = 300)

####################
#Make mutation frequency over time plots
####################
#make plot function
plot_temp_dyn_of_adap_change <- function(variant_table_dat, cand_adap_site_info, plot_title){
	variant_table_dat %>% mutate(passage = recode_factor(passage, "Clinical sample" = "P0")) %>%
	ggplot(aes(x = passage, y = fMut)) +
	geom_point(aes(col = variant_sample, alpha = alpha), size = 1) + 
	geom_line(aes(group = exp_id_pos, col = variant_sample, linetype = replicate, alpha = alpha), size = 0.5) + 
	
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
		name = "Replicate",
		values = c( 
			"A" = "solid", 
			"B" = "dashed"
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

	geom_hline(yintercept = fMut_lw_threshold, linetype = "dashed", size = 0.5) +

	scale_alpha(guide = 'none') + 

	scale_x_discrete(name = "Passage") +
	ylim(0,1) + 
	scale_y_continuous(name = "Mutation frequency") +
	facet_wrap(
		. ~ position, 
		ncol = 6, 
		labeller = as_labeller(
			setNames(
				cand_adap_site_info$strip_lab,
				cand_adap_site_info$position
			)
		)
	) +
	labs(title = plot_title) +

#	geom_text(
#		data = cand_adap_site_info %>% 
#			mutate(
#				label_B.1.36.16 = 
#				ifelse(
#					adap_site_cand_B.1.36.16, 
#					paste(
#						"B.1.36.16: \u0394s = ", round(mut_adv_coef_B.1.36.16, 3), 
#						";\n\u00B5 = ", mu_nt_aa_B.1.36.16, 
#						";\n+ve experiment = ", positive_exp_id_B.1.36.16,
#						sep = ""
#					),
#					NA
#				),
#				label_AY.30 = 
#				ifelse(
#					adap_site_cand_AY.30, 
#					paste(
#						"AY.30: \u0394s = ", round(mut_adv_coef_AY.30, 3), 
#						";\n\u00B5 = ", mu_nt_aa_AY.30, 
#						";\n+ve experiment = ", positive_exp_id_AY.30,
#						sep = ""
#					),
#					NA
#				)
#			) %>% 
#			unite(label, label_B.1.36.16, label_AY.30, na.rm = TRUE, sep = "\n", remove = FALSE) %>% 
#			select(position, label),
#		mapping = aes(x = 0, y = 1, label = label),
#		size = 1.75,
#		hjust = 0, vjust = 1, 
#	) + 

	theme_classic() +
	theme(
		strip.text = element_text(size = 6),
		axis.text = element_text(size = 6),
		axis.title = element_text(size = 8, face = "bold"),
		
		plot.title = element_text(size = 10, face = "bold"),
	)
}

#make Ext Data Figure S3-S6: Temporal dynamics of potential adaptive changes found in all cell lines
mu_freq_plots <- NULL
mu_freq_plots[["Vero E6"]] <- plot_temp_dyn_of_adap_change(
	variant_table_dat = variant_table_dat_for_plotting %>% 
		filter(cell_line == "Vero E6") %>%
		mutate(variant_sample = paste(variant, sample, sep = ": ")),
	cand_adap_site_info = cand_adap_site_info %>%
		filter(cell_line == "Vero E6") %>%
		mutate(subdomain = ifelse(subdomain == "", NA, subdomain)) %>%
		unite(strip_lab, gene, subdomain, na.rm = TRUE, sep = "\n", remove = FALSE) %>%
		mutate(strip_lab = paste(prettyNum(position, big.mark = ","), strip_lab, sep = ": ")), 
	plot_title = "Temporal dynamics of potential Vero E6 adaptive mutations") + 
	theme(
		legend.position = "bottom",
		legend.box = "horizontal",
		legend.key.width = unit(0.75,"cm"), 
		legend.key.height = unit(0.25,"cm"), 
		legend.margin = margin(t = 0, unit = "cm"),
		legend.title = element_text(size = 8, face = "bold"), 
		legend.text = element_text(size = 6),
	) + 
	guides(colour = guide_legend(title.position = "top"),
      	linetype = guide_legend(title.position = "top"))

mu_freq_plots[["Vero E6-TMPRSS2"]] <- plot_temp_dyn_of_adap_change(
	variant_table_dat = variant_table_dat_for_plotting %>% 
		filter(cell_line == "Vero E6-TMPRSS2") %>%
		mutate(variant_sample = paste(variant, sample, sep = ": ")),
	cand_adap_site_info = cand_adap_site_info %>%
		filter(cell_line == "Vero E6-TMPRSS2") %>%
		mutate(subdomain = ifelse(subdomain == "", NA, subdomain)) %>%
		unite(strip_lab, gene, subdomain, na.rm = TRUE, sep = "\n", remove = FALSE) %>%
		mutate(strip_lab = paste(prettyNum(position, big.mark = ","), strip_lab, sep = ": ")), 
	plot_title = "Temporal dynamics of potential Vero E6/TMPRSS2 adaptive mutations") + 
	theme(
		legend.position = "bottom",
		legend.box = "horizontal",
		legend.key.width = unit(0.75,"cm"), 
		legend.key.height = unit(0.25,"cm"), 
		legend.margin = margin(t = 0, unit = "cm"),
		legend.title = element_text(size = 8, face = "bold"), 
		legend.text = element_text(size = 6),
	) + 
	guides(colour = guide_legend(title.position = "top"),
      	linetype = guide_legend(title.position = "top"))


mu_freq_plots[["Calu-3"]] <- plot_temp_dyn_of_adap_change(
	variant_table_dat = variant_table_dat_for_plotting %>% 
		filter(cell_line == "Calu-3") %>%
		mutate(variant_sample = paste(variant, sample, sep = ": ")),
	cand_adap_site_info = cand_adap_site_info %>%
		filter(cell_line == "Calu-3") %>%
		mutate(subdomain = ifelse(subdomain == "", NA, subdomain)) %>%
		unite(strip_lab, gene, subdomain, na.rm = TRUE, sep = "\n", remove = FALSE) %>%
		mutate(strip_lab = paste(prettyNum(position, big.mark = ","), strip_lab, sep = ": ")), 
	plot_title = "Temporal dynamics of potential Calu-3 adaptive mutations") + 
	theme(
		legend.position = "bottom",
		legend.box = "horizontal",
		legend.key.width = unit(0.75,"cm"), 
		legend.key.height = unit(0.25,"cm"), 
		legend.margin = margin(t = 0, unit = "cm"),
		legend.title = element_text(size = 8, face = "bold"), 
		legend.text = element_text(size = 6),
	) + 
	guides(colour = guide_legend(title.position = "top"),
      	linetype = guide_legend(title.position = "top"))


mu_freq_plots <- align_plots(plotlist = mu_freq_plots, align = "hv", axis = "tblr")

#save figure to files
ggsave(Ext_Dat_FigS3_temp_dyn_of_adap_changes.E6.png_filename,
	plot = mu_freq_plots[["Vero E6"]],
	width = 16, height = 63.33, units = "cm",
	dpi = 300)

ggsave(Ext_Dat_FigS3_temp_dyn_of_adap_changes.E6.svg_filename,
	plot = mu_freq_plots[["Vero E6"]],
	width = 16, height = 63.33, units = "cm",
	dpi = 300)

ggsave(Ext_Dat_FigS4_temp_dyn_of_adap_changes.TM.png_filename,
	plot = mu_freq_plots[["Vero E6-TMPRSS2"]],
	width = 16, height = 9.80, units = "cm",
	dpi = 300)

ggsave(Ext_Dat_FigS4_temp_dyn_of_adap_changes.TM.svg_filename,
	plot = mu_freq_plots[["Vero E6-TMPRSS2"]],
	width = 16, height = 9.80, units = "cm",
	dpi = 300)

ggsave(Ext_Dat_FigS5_temp_dyn_of_adap_changes.C3.png_filename,
	plot = mu_freq_plots[["Calu-3"]],
	width = 16, height = 46.86, units = "cm",
	dpi = 300)

ggsave(Ext_Dat_FigS5_temp_dyn_of_adap_changes.C3.svg_filename,
	plot = mu_freq_plots[["Calu-3"]],
	width = 16, height = 46.86, units = "cm",
	dpi = 300)

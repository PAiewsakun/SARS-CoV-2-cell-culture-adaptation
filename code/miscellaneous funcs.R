#make ref genome structure plot
make_ref_genome_structure_plot <- function(ref_annotations_filename = "data/ref genome annotation.txt", text_size = 2.5, h = 0.05){
	#load curated ref SARS-CoV-2 Wuhun Hu 1 genome annotation 
	ref_annotations <- read.table(ref_annotations_filename, header = TRUE, sep = "\t", quote = "")
	genome_length <- max(ref_annotations$end)
	ref_annotations <- ref_annotations %>% filter(type == "gene") %>% mutate(start = as.numeric(start), end = as.numeric(end)) 

	#process gene annotation data
	#compute gene-box mid locations
	ref_annotations <- ref_annotations %>% mutate(mid = rowMeans(cbind(start, end)), .after = end)

	#add gene type info to the data
	ref_annotations <- ref_annotations %>% 
		mutate(
			type = ifelse(
				grepl("NSP", gene, fixed = TRUE), 
				"NSP coding gene", 
				ifelse(
					grepl("ORF", gene, fixed = TRUE),
					"Accessory gene",
					"Structural gene"
				)
			) 
		) %>% 
		arrange(type)

	#set gene box colours
	ref_annotations <- ref_annotations %>% 
		mutate (
			fill = c(
				rep_len(c("gold", "gold3"),sum(type=="Accessory gene")),
				rep_len(c("steelblue3", "steelblue1"),sum(type=="NSP coding gene")),
				rep_len(c("firebrick3", "firebrick1"),sum(type=="Structural gene"))
			)
		) %>% 
		arrange(start)

	#make the plot
	ref_genomic_str_plot <- ggplot() + 
		geom_segment(aes(x = 1, y = 0, xend = genome_length, yend = 0)) +
		geom_tile(
			data = ref_annotations,
			aes(	
				x = mid, 
				y = 0, 
				width = abs(start - end), 
				height = h
			), 
			fill = ref_annotations[,"fill"],
			colour = NA
		) +
		geom_text(
			data = ref_annotations,
			aes(x = mid, y = 0, label = gene,),
			col = ref_annotations[,"fill"],
			hjust = 0, nudge_y = (h/2+0.01), 
			angle = 90, 
			size = text_size, fontface = 3
		) +
		scale_x_continuous(limits = c(1,29903), expand = c(0, 0)) +
		scale_y_continuous(, limits = c(-h, h + 0.1)) + 
		theme_void() 
	return(ref_genomic_str_plot)
}

#data pre processing functions
set_factors <- function(dat){
	dat <- dat %>% mutate (	
		variant = factor(variant, levels = c("B.1.36.16", "AY.30")),
		sample = factor(sample, levels = c("73NLt", "CV130", "OTV54", "NH783")),
		cell_line = factor(cell_line, levels = c("Clinical sample", "Vero E6", "Vero E6-TMPRSS2", "Calu-3")),
		passage = factor(passage, levels = c("Clinical sample", "P1", "P2", "P3", "P4"), ordered = TRUE),
		replicate = factor(replicate, levels = c("A", "B")),

		variant_sample = factor(paste(variant, sample, sep = ": "), levels = c("B.1.36.16: 73NLt", "B.1.36.16: CV130", "AY.30: OTV54", "AY.30: NH783")) #add variant_sample column
	) 

	dat <- dat %>% mutate_if(is.character, as.factor)
	return(dat)
}

#get gene names from position
get_gene_from_pos <- function(pos, ref_annotations){
	ref_annotations$gene[unlist(lapply(sapply(pos, function(x) which(ref_annotations$start <= x & x <= ref_annotations$end)), function(x) ifelse(length(x)==0, NA, x) ))]
}

#get gene subregion names from position
get_region_from_pos <- function(pos, ref_annotations){
	ref_annotations$region[unlist(lapply(sapply(pos, function(x) which(ref_annotations$start <= x & x <= ref_annotations$end)), function(x) ifelse(length(x)==0, NA, x) ))]
}

#remove empty strings from a vector
remove_empty_strings_from_a_vector <- function(x){
	x[x != ""]
}

#deduce aa mutation from nt mutation
nucleotide_mu_to_protein_mu <- function(position, ori, mu, ref, ref_annotations){
	start <- ref_annotations %>% filter(start <= position, position <= end) %>% pull(start)
	end <- ref_annotations %>% filter(start <= position, position <= end) %>% pull(end)
	aa_pos <- floor((position - start)/3)+1

	ori_seq <- ref
	ori_seq[position] <- ori
	ori_aa <- as.character(trans(as.DNAbin(ori_seq[start:end]))[aa_pos])# ape::trans; ape::as.DNAbin

	mu_seq <- ref
	mu_aa <- NULL
	for (m in mu){
		m = ifelse(m == "DEL", "*", ifelse(m == "INS", "*" , m))
		mu_seq[position] <- m
		mu_aa <- c(mu_aa, as.character(trans(as.DNAbin(mu_seq[start:end]))[aa_pos]))
	}

	mu_aa[mu_aa == ori_aa] <- ""
	return(paste(ori_aa, prettyNum(aa_pos, big.mark = ","), toupper(mu_aa), sep = ""))
}

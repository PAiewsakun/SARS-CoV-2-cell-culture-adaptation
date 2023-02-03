# SARS-CoV-2-cell-culture-adaptation

03/02/2023

Description
==============
R scripts used to generate figures and results in the manuscript entitled "Systematic exploration of SARS-CoV-2 adaptation to Vero E6, Vero E6/TMPRSS2, and Calu-3 cells" (Aiewsakun et al. 2022).

The analyses were run using the following R packages:
==============
ape == 5.6.2,
cowplot == 1.1.1,
dplyr == 1.0.9,
forcats == 0.5.1,
ggh4x == 0.2.1,
ggnewscale == 0.4.5,
ggplot2 == 3.3.6,
hash == 2.2.6.1,
lme4qtl == 0.2.2,
R.utils == 2.11.0,
stringr == 1.4.0,
tibble == 3.1.7,
tidyr == 1.2.0,
tidyverse == 1.3.1,

in R 4.0.4

INSTALLATION INSTRUCTIONS
==============
Simply download the repo. If necessary, install required packages in R with: install_packages(c("ape", "cowplot", "dplyr", "forcats", "ggh4x", "ggnewscale", "ggplot2", "hash", "R.utils", "stringr", "tibble", "tidyr", "tidyverse")). To install "lme4qtl", see https://github.com/variani/lme4qtl

INSTRUCTIONS FOR USE
==============
Upzip the file "variant table.zip" in the "data" directory, and leave the file "variant table.txt" there.
Run the scrtips in R.
To execute the codes, set the variable "path_to_wd" to the path of the "code" directory.
Results will be available to you in "code/results".

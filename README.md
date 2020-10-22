This repository contains experimental software, data and analysis code for the results reported in:

Popov, V., So, M. & Reder, L. (2020). Word frequency affects binding probability not memory precision. PsyArXiv.

The code runs with paths relative to the parent folder, which is set-up in every script by the here() package and the setwd(here()) command, which finds the .Rproj file and sets that as the working dir.

All code runs under the following R, OS and packages versions:

```
R version 4.0.2 (2020-06-22)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19041)

Matrix products: default

Random number generation:
 RNG:     Mersenne-Twister 
 Normal:  Inversion 
 Sample:  Rounding 
 
locale:
[1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252    LC_MONETARY=English_United States.1252
[4] LC_NUMERIC=C                           LC_TIME=English_United States.1252    

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] lme4_1.1-23     Matrix_1.2-18   circular_0.4-93 forcats_0.5.0   stringr_1.4.0   dplyr_1.0.2     purrr_0.3.4     readr_1.3.1     tidyr_1.1.2    
[10] tibble_3.0.3    ggplot2_3.3.2   tidyverse_1.3.0
```
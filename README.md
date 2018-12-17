# qtl2pleio-manuscript
Resources from the qtl2pleio manuscript, "Testing pleiotropy vs. separate QTL in multiparental populations"


There are 4 subdirectories in this repository:

1. Rmd  
2. chtc  
3. R  
4. overleaf-repo

We include a Makefile to ease manuscript assembly. Just type `make` at the BASH prompt to compile the pdf from the tex file "overleaf-repo/main.tex".




## R contents

Subdirectory "R" contains R scripts and outputted figures.




## Rmd contents

Rmd contains Rmarkdown files and github markdown documents (ie, the output from the Rmarkdown documents). It also contains saved png files for use in the manuscript. Additionally, intermediate results files are here, namely, those that contain the `pvl_scan` results for each pair in the Recla analysis. 

## overleaf-repo contents

overleaf-repo contains files that we originally created on the website overleaf.com. They include the tex file "overleaf-repo/main.tex" for compilation with latex.


## chtc contents

Much of the computing was done through the Center for High-Throughput Computing at UW-Madison. I have placed the needed files for use with the chtc in this subdirectory. 

There are 3 subdirectories:

1. Recla-bootstrap
2. type-I-error-rate-study
3. power-study


## Contents of chtc/Recla-bootstrap

Within Recla-bootstrap, there are 5 subdirectories: 

1. submit_files  
2. Rscript  
3. shell_scripts  
4. squid
5. data

Each contains 3 files, all of which we used in analysis of the Recla data. Additional files are needed to fully reproduce our analysis. For example, several files are placed on a squid server for fast transfer when using the CHTC. We've placed the needed files in the "squid" subdirectory. We've also placed the necessary data files in the subdirectory "data".


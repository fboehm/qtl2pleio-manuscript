RMD_FILE = Rmd/2018-11-17_recla-figures-revised.Rmd

# Tools
LATEXMK = latexmk
RM = rm -f

# Project-specific settings
DOCNAME = overleaf-repo/main

# Targets
all: doc
doc: pdf
pdf: overleaf-repo/main.pdf
fig: R/power-curves.eps Rmd/scatter.eps Rmd/chr8-lods.eps Rmd/coefs.eps Rmd/profile.eps Rmd/genomewide_lod_trait10.eps Rmd/genomewide_lod_trait22.eps


# Rules
%.pdf: %.tex fig
	cd overleaf-repo; $(LATEXMK) -pdfps -bibtex main.tex
mostlyclean:
	cd overleaf-repo; $(LATEXMK) -silent -c
	$(RM) overleaf-repo/*.bbl

clean: mostlyclean
	cd overleaf-repo; $(LATEXMK) -silent -C
	$(RM) overleaf-repo/*.run.xml overleaf-repo/*.synctex.gz
	$(RM) overleaf-repo/*.d

.PHONY: all clean doc mostlyclean pdf

# Include auto-generated dependencies
#-include *.d


R/power-curves.eps: R/plot-power-curves.R
	R CMD BATCH R/plot-power-curves.R
	
Rmd/scatter.eps: $(RMD_FILE)
	Rscript -e "rmarkdown::render('$<')"
  

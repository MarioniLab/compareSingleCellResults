targets: all

RCMD=Rdevel

# Obtaining the vignette sources.
src:
	 git clone https://github.com/MarioniLab/compareSingleCell src

update: src
	cd src && git pull
	${RCMD} CMD INSTALL --preclean src/

# Moving the scripts to avoid directory chaos.
%.Rmd: src
	cp src/vignettes/$@ .

ref.bib: src
	cp src/vignettes/ref.bib .

# Creating the *.knit.md file (destroying it if the command fails).
%.knit.md : %.Rmd ref.bib
	${RCMD} --no-save --slave -e "rmarkdown::render('$<', clean=FALSE)" || rm $@

embryo_merge.knit.md: embryo_preprocess.knit.md

embryo_abundance.knit.md: embryo_merge.knit.md

embryo_expression.knit.md: embryo_merge.knit.md

all: embryo_preprocess.knit.md \
	embryo_merge.knit.md \
	embryo_abundance.knit.md \
	embryo_expression.knit.md 

# Cleaning commands.
clean: 
	rm -rf *.Rmd *.html *_files *_cache *.knit.md 

distclean: clean
	rm -rf raw_data
	rm -rf src

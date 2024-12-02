
######################
######################
# makefile generated 
# dude makefile 
# created by s.e.pav 
# $Id: Makefile 89 2006-01-26 20:59:08Z spav $
######################
######################

############### FLAGS ###############

# these don't work at times
# in that case set these by hand?

LATEX       := $(shell which latex)
BIBTEX      := $(shell which bibtex)
PDFLATEX    := $(shell which pdflatex)
MAKEINDEX   := $(shell which makeindex)
RSCRIPT 		:= $(shell which Rscript)
LITTLER 		:= $(shell which littler)

RLIB         = /usr/lib64/R

TEXINPADD    = .:./Definitions:$(RLIB)/share/texmf/tex/latex

PRETEX       = TEXINPUTS=$(TEXINPADD):$$TEXINPUTS
PREBIB       = BSTINPUTS=$(TEXINPADD):$$BSTINPUTS \
               BIBINPUTS=$(TEXINPADD):$$BIBINPUTS 

PREIDX       = INDEXSTYLE=$(TEXINPADD):$$INDEXSTYLE

PROJECT      = 

#SAVE
# tracked projects
PROJECTS     = fromo
ALL_PDF 		 = $(foreach proj,$(PROJECTS),$(proj).pdf)

#UNSAVE
# add on dependencies (subchapters of snrinf)
R_DEPS 			 = 
TEX_EXTRAS   = sep.sty common.bib
# nonlocal dependencies
STY_FILES    = 

#aspell
ASPELL_FLAGS = 

# for git
VERSION 					?= 2

DO_TIME 			?= 

ifeq ($(DO_TIME),1)
TIMEIT 				= time
else
TIMEIT 				= 
endif

############## DEFAULT ##############

default : all

############## MARKERS ##############

.PHONY   : 
.SUFFIXES: .tex .bib .dvi .ps .pdf .eps
.PRECIOUS: %.pdf 

############ BUILD RULES ############

.PHONY   : help 

# this will have to change b/c of inclusion file names...
help:  ## generate this help message
	@grep -h -P '^(([^\s]+\s+)*([^\s]+))\s*:.*?##\s*.*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}'


.PHONY   : all

# an easy target
all : $(ALL_PDF) ## build the document by knitting source code

%.tex : %.Rnw $(R_DEPS)
		$(TIMEIT) $(RSCRIPT) -e 'require(knitr);knit("$<")'

%.R : %.Rnw
		$(TIMEIT) $(RSCRIPT) -e 'require(knitr);knit("$<",tangle=TRUE)'

%.pdf : %.qmd
		$(TIMEIT) $(RSCRIPT) -e 'require(quarto);quarto_render("$<",output_format="beamer")'

%.pdf : %.tex
	latexmk -f -bibtex -pdf -pdflatex="$(PDFLATEX)" -use-make $<

# tex extras
%.bbl : %.bib
	$(PREBIB) $(BIBTEX) $*

%.bbl : %.aux
	$(PREBIB) $(BIBTEX) $*

%.ind : %.idx
	$(PREIDX) $(MAKEINDEX) $*

############# CLEAN UP ##############

# clean up
%.clean : 
		-rm -f $*.aux $*.log $*.dvi $*.bbl $*.blg $*.toc $*.ilg $*.ind
		-rm -f $*.out $*.idx $*.lot $*.lof $*.brf $*.nav $*.snm
%.realclean : %.clean
		-rm -f $*.ps $*.pdf
		-rm -f $*-[0-9][0-9][0-9]*.eps $*-[0-9][0-9][0-9]*.pdf

############### RULES ###############


release.tex: snrinf.tex
	perl -pe 's{figure/}{};' < $< > $@

release : release.tex
	mv release.tex snrinf.tex
	@-echo "upload to arxiv"

# check it

# clean up
clean: $(patsubst %,%.clean,$(PROJECTS))
realclean: $(patsubst %,%.realclean,$(PROJECTS))
		-rm -f Rplots.pdf

cleancache: 
		echo "killing knitr cache! ack!"
		-rm -rf cache

superclean: realclean cleancache

.tags :
	nice -n 18 ctags -f .tmp_tags --recurse --language-force=R --fields=+i `find . -regextype posix-egrep -regex '.*.R(nw)?'`;
	mv .tmp_tags $@

######################
######################

tag : ## advice on github tagging
	@-echo "git tag -a v$(VERSION) -m 'release $(VERSION)'"
	@-echo "git push --tags"

untag : ## advice on github untagging
	@-echo "git tag --delete v$(VERSION)"
	@-echo "git push origin :v$(VERSION)"

really_tag : ## actually github tag 
	git tag -a v$(VERSION) -m 'release $(VERSION)'
	git push --tags

#for vim modeline: (do not edit)
# vim:ts=2:sw=2:tw=149:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:tags=tags;:syn=make:ft=make:ai:si:cin:nu:fo=croqt:cino=p0t0c5(0:

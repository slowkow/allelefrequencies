.PHONY: all

all: README.md afnd.tsv

README.md: README.Rmd
	R -e 'rmarkdown::render("README.Rmd")'

afnd.tsv: allelefrequencies.py
	python allelefrequencies.py

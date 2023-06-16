.PHONY: all

all: README.md README.html afnd.tsv

README.md: README.Rmd
	R -e 'rmarkdown::render("README.Rmd", "all")'

README.html: README.Rmd
	R -e 'rmarkdown::render("README.Rmd", "all")'

afnd.tsv: allelefrequencies.py
	python allelefrequencies.py

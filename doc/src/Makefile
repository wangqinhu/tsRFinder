UNAME := $(shell uname)

pdf: manual.tex
	pdflatex ./manual.tex
	pdflatex ./manual.tex
	mv manual.pdf ..
	rm ./manual.out ./manual.aux ./manual.toc ./manual.log

ifeq ($(UNAME), Darwin)
	open ../manual.pdf
endif

clean:
	rm ../manual.pdf

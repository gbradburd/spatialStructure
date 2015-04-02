%.html : %.md
	pandoc $< > $@

%.html : %.Rmd
	R --vanilla -e "require(rmarkdown);render(\"$<\",output_file=\"$@\")"

%.png : %.pdf
	convert -density 300 $< -flatten $@

%.svg : %.ink.svg
	inkscape $< --export-plain-svg=$@

%.pdf : %.ink.svg
	inkscape $< --export-pdf=$@

%.svg : %.pdf
	inkscape $< --export-plain-svg=$@

%.png : %.pdf
	convert -density 300 $< -flatten $@


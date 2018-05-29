all:  
	#xsltproc  --stringparam html.stylesheet docbook.css mystyles.xsl index.xml
	xsltproc  --stringparam html.stylesheet style.css mystyles.xsl index.xml
	#xsltproc --stringparam html.stylesheet docbook.css mystyles.xsl index.xml
	runsed *.html
	mv *.html ../html
	cp *.css ../html
	#
	# This creates the README file that goes in the top level directory 
	docbook2txt index.xml
	mv index.txt ../../README

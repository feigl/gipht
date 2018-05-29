<?xml version='1.0'?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
                version='1.0'
                xmlns="http://www.w3.org/TR/xhtml1/transitional"
                exclude-result-prefixes="#default">

<!-- this import is probably machine and/or distribution dependent. -->
<!-- Change to "docbook.xsl" in order to get an all-in-one xhtml page. -->
<!-- 
<xsl:import href="/usr/share/sgml/docbook/xsl-stylesheets-1.64.0/xhtml/docbook.xsl"/>
<xsl:import href="/opt/matlab7/sys/namespace/docbook/v4/xsl/html/docbook.xsl"/>
<xsl:import href="/usr/share/apps/ksgmltools2/docbook/xsl/html/docbook.xsl"/>
-->
<xsl:import href="/usr/share/sgml/docbook/xsl-stylesheets-1.69.1/xhtml/chunk.xsl"/>

<!-- Set the stylesheet -->
<xsl:param name="html.stylesheet" select="'docbook.css'"/>

<xsl:param name="toc.section.depth" select="1"/>
<xsl:param name="generate.section.toc.level" select="1"/>
<xsl:param name="chunk.section.depth" select="2"/>
<!-- 
<xsl:param name="section.autolabel" select="1"/>
<xsl:param name="xref.with.number.and.title" select="1"/>
Set the stylesheet -->



</xsl:stylesheet>

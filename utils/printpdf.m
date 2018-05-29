function printpdf(pdfFileName)

% Save the pdf (this is the same method used by "saveas")
dpi = 600;
print(gcf,'-dpdf',pdfFileName,sprintf('-r%d',dpi));
return

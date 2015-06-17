[Gx,Gy,Gv]= grdread2('Quad_ASC_CSK_2013.75_2014_04_09_lonlatvlos.grd');
figure;
imagesc(Gx,Gy,Gv*1000);axis xy
cmapblacknan
h=colorbar('Location','SouthOutside');
xlabel(h,'mm/yr');
title(h,'LOS velocity');
printpdf('Quad_ASC_CSK_2013.75_2014_04_09_lonlatvlos.pdf');

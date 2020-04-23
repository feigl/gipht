% test printpdf function
figure;
bar([1 10 7 8 2 2 9 3 6])
title('title string')
xlabel('x label');
ylabel('y label');
printpdf('test1.pdf',gcf,600)
printpdf('test2.pdf',gcf,600,true)
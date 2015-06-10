function bpest = adjustbp(tepochs,DD,bp,species, iuniqorbs, uniqdates)
%
% given time tags and Bperp values, perform temporal adjustment
%
% bpest = adjustbp(tepochs,DD,bp,species)
%  
% bpest = adjustbp(tepochs,DD,bp,species, iuniqorbs, uniqdates)
%
% Kurt Feigl CNRS 
% 2005 January
% 2005 JUL 14 add orbit numbers and dates as option
%
% 2007 DEC replace Bp with Ddop


fprintf(1,'%s begins ...\n',mfilename);

nargchk(4,6,nargin);
nargoutchk (1,1,nargout);


if (max(bp)-min(bp) < 1)
   ylab = 'Dopp/PRF';
else
   ylab = 'Bperp(m)';
end

% design matrix
DD2 = DD;

% data vector
data = bp;

% number of data
nd = length(bp);

% number of species
nm = size(species); nspecies = nm(1);

% number of pairs
nm = size(DD); np = nm(1);

% number of epochs
me = nm(2);

% add constraining equations as lines to DD2
for j=1:nspecies
    family = species(j,:);
    k=isfinite(family);
    mk = length(find(k == 1));
    family = family(1:mk);
    nd = nd+1;
    for i = 1:me
        DD2(nd,i) = 0;
    end
    for i = 1:mk
        DD2(nd,family(i)) = 1;
        data(nd) = 0;
    end
end


disp 'Length of data vector, including constraints';md = size(data)
disp 'Dimensons of design matrix, including constraints';nmdd2 = size(DD2)
disp 'Rank defiency, including constraints'; rd = me - rank(DD2)

disp 'Begin least squares adjustment...'

if rd <= 1 
   if nd > 2
      bpest = inv(DD2' * DD2) * (DD2' * data);
   else
      bpest(1) = 0;  % arbitrarily make first the origin
      bpest(2) = data(1);
   end
else
   if nd > 2
      warning('Rank defiency persists. Using pseudoinverse!');
      bpest = pinv(DD2' * DD2) * (DD2' * data);
   else
      bpest(1) = 0;  % arbitrarily make first the origin
      bpest(2) = data(1);
   end
end

fid=fopen(sprintf('%s.bp',mfilename),'w');
if nargin < 6
	fprintf(1  , 'Iepoch  Year   %s\n',ylab);
	fprintf(fid, 'Iepoch  Year   %s\n',ylab);
    
	for i=1:me
           fprintf(1,  '%3d %11.4f %20.4f\n',i,tepochs(i),bpest(i));
           fprintf(fid,'%3d %11.4f %20.4f\n',i,tepochs(i),bpest(i));
	end
else
	% tell us about the individual epochs 
	fprintf(1,   'index orbnum date     species yr %s\n',ylab);
	fprintf(fid, 'index orbnum date     species yr %s\n',ylab);
	for i = 1:me
           [j,k] = find(species == i);
           j=j(1);k=k(1);
          %fprintf (1,  '%5d %c %7d %s %3d %10.4f %10.4f\n',i,char(i+64),iuniqorbs(i),uniqdates{i},j,tepochs(i),bpest(i));
          %fprintf (fid,'%5d %c %7d %s %3d %10.4f %10.4f\n',i,char(i+64),iuniqorbs(i),uniqdates{i},j,tepochs(i),bpest(i));
           fprintf (1,  '%5d %7d %s %c %10.4f %10.4f\n',i,iuniqorbs(i),uniqdates{i},char(j+64),tepochs(i),bpest(i));
           fprintf (fid,'%5d %7d %s %c %10.4f %10.4f\n',i,iuniqorbs(i),uniqdates{i},char(j+64),tepochs(i),bpest(i));
	end
end
fclose(fid);



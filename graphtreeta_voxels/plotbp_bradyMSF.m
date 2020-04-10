function ktours = plotbp_sanemMSF(tepochs, bpest, DD, species, iuniqorbs, uniqdates, plotts,ylab, sat_id, cal_date)
%function ktours = plotbp(tepochs, bpest, DD, trees, iuniqorbs, uniqdates, plotts,ylab)
%
% plot pseudo-absolute Bperp as a function of time 
% and return ktours the mininum-lenthg traveling salesman path
%
% Automatically handle 3 cases
% 1  X = t      Y = Bperp
% 2  X = t      Y = Ddop
% 3  X = Bperp  Y = Ddop
%
%
% ktours = plotbp(tepochs, bpest, DD, species, iuniqorbs, uniqdates, plotts)
%
%
% ktours = plotbp(tepochs, bpest, DD, species)
%  
% ktours = plotbp(tepochs, bpest, DD, species, iuniqorbs, uniqdates)
%
% ktours = plotbp(tepochs, bpest, DD, species, iuniqorbs, uniqdates, plotts)
%     plotts = 0 do not connect points with traveling salesman net 
%     plotts = 1 connect points with traveling salesman net 
%     plotts = 2 connect points with traveling salesman net AND existence
%     plotts = 3 connect points with existence only
%
% Kurt Feigl CNRS 
% 2005 January
% 2005 JUL 14 add orbit numbers and dates as option
% 2006 FEB 23 fix mod bug
% 2007 NOV 17 connect points with traveling salesman trajectory
% 2008-MAR-29 correct annoying bug
% 2014-JUL-06 include ylab as input
%
% Elena C. Baluyut, UW-Madison
% 2015-01-30 add text to plot showing month and day of epochs per interval

FigHandle = figure;

fidtxtout = fopen(sprintf('%sout.txt',mfilename),'a+t');
for ifile = [1 fidtxtout]
   fprintf(ifile,'%s begins at %s\n',mfilename,datestr(now,31));
end

%set(FigHandle,'defaulttextinterpreter','latex');
nargchk(4,8,nargin);
if nargin <= 6
   plotts = 0;
end

if nargin == 10
    lab_on = 1;
else
    lab_on = 0;
end

nargoutchk (1,1,nargout);

if (max(bpest)-min(bpest) < 1.0)
   if max(tepochs) < 1990
      icase = 3;
      
      xlab = '$B_{\perp}$ (m)';
      ylab = 'Azimuthal Doppler/PRF';
      Tend = 10;
      yrbp = 100/0.1; % 100 m of Bperp is like 0.1 PRF of Doppler separatio
   else
      icase = 2;
      xlab = 'year';
      ylab = 'Azimuthal Doppler/PRF';
      Tend  = 30;
      yrbp = 1/0.001; % 1 year of time separation is like 0.001 PRF of Doppler separation
   end
else
    if exist(ylab,'var') == 1
        icase = 1;
        xlab = 'year';
        Tend = 25;
        yrbp = 1; % 1 year of time separation is like 1000 m of orbits separation
    else
        icase = 1;
        xlab = 'year';
        ylab = '$B_{\perp}$ (m)';
        Tend = 25;
        yrbp = 1/1000; % 1 year of time separation is like 1000 m of orbits separation
    end
end

if plotts == 1
   titl = 'Optimal set of pairs';
else
   titl = 'Selected pairs';
end

% define symbols to use

%            y     yellow        .     point              -     solid
%            m     magenta       o     circle             :     dotted
%            c     cyan          x     x-mark             -.    dashdot 
%            r     red           +     plus               --    dashed   
%            g     green         *     star
%            b     blue          s     square
%            w     white         d     diamond
%            k     black         v     triangle (down)
%                                ^     triangle (up)
%                                <     triangle (left)
%                                >     triangle (right)
%                                p     pentagram
%                                h     hexagram

% mysyms = {'gx-' 'ro-' 'b*-' 'ks-' 'md-' 'cv-'};
% mysols = {'g-'  'r-'  'b-'  'k-'  'm-'  'c-'};
% mylins = {'g:'  'r:'  'b:'  'k:'  'm:'  'c:'};
% mydash = {'g--' 'r--' 'b--' 'k--' 'm--' 'c--'};
% mysym0 = {'gx'  'ro'  'b*'  'ks'  'md'  'cv'};
mysyms = {'gx-' 'ro-' 'b*-'  };
mysols = {'g-' 'r-'  'b-'  };
mylins = {'g:' 'r:'  'b:' };
mydash = {'g--' 'r--' 'b--'};
mysym0 = {'gx' 'ro'  'b*'};



% graphics handle to return
%h=figure; 

% number of species
nm = size(species); nf = nm(1);

% number of pairs
nm = size(DD); np = nm(1);

%plot origin to make legend come out right
for j=1:nf
    plot(min(tepochs),0,mysyms{1+mod(j,length(mysyms))});
    hold on;
    specie = species(j,:);
    k=isfinite(specie);
    specie=specie(k);
    me = length(find(k == 1));
    if nargin >= 6
       famnam{j} = strcat(sprintf('Component %s epochs:',char(j+64)),sprintf(' %1d',iuniqorbs(specie(1:me))));
    else
       famnam{j} = strcat(sprintf('Component %s ID:',char(j+64)),sprintf('%3d',char(j+64),specie(1:me)));
    end
end
ktours = zeros(size(species));

% legend(famnam,'Location','NorthOutside');

% over plot origin with white
plot(min(tepochs),0,'sw');
plot(min(tepochs),0,'ow');
plot(min(tepochs),0,'xw');
plot(min(tepochs),0,'*w');
plot(min(tepochs),0,'dw');
plot(min(tepochs),0,'vw');

% draw end points of available pairs
for i=1:np      
   ddcol = DD(i,:); 
   j=find(abs(ddcol)>0);
   id0(i) = min(j);  
   id1(i) = max(j);
end
id2 = unique([id0 id1]);


% for ifile = [1 fidtxtout]
%    fprintf(ifile,'Pair Component Member0 Member1 orbn0 orbn1 year0 year1 %s %s\n',xlab,ylab);
% end
i=0;
nt = numel(tepochs);

if plotts > 0
   for j = 1:nf
      % find the members of this species
      kkeep = find(isnan(species(j,:)) == 0);
      tspxy = zeros(numel(kkeep),2);
      % traveling salesman coordinates are time and Bperp
      tspxy(:,1) = tepochs(species(j,kkeep));
      tspxy(:,2) = bpest(species(j,kkeep));

      % rescale
      tspxy(:,1) = tspxy(:,1)*yrbp;
      nmem = numel(kkeep);
      
      if nmem > 3
         % traveling salesman problem
         ktour = tspsiman(tspxy,Tend);
         ktour = ktour(1:length(ktour)-1); % remove -1
      else
         if nmem == 3
            dista = sqrt( (tspxy(2,1))-(tspxy(1,1))^2 + (tspxy(2,2))-(tspxy(1,2))^2);
            distb = sqrt( (tspxy(3,1))-(tspxy(2,1))^2 + (tspxy(3,2))-(tspxy(2,2))^2);
            distc = sqrt( (tspxy(3,1))-(tspxy(1,1))^2 + (tspxy(3,2))-(tspxy(1,2))^2);
            if dista + distb < dista + distc
               ktour = [1 2 3];
            else
               ktour = [1 3 2];
            end
         else
            ktour = [1 2];
         end
      end
      ktours(j,1:numel(ktour)) = ktour;
 
      for k=1:numel(ktour)-1
         i=i+1;  % count pairs
         i0=species(j,kkeep(ktour((k))));
         i1=species(j,kkeep(ktour((k+1))));
         % Traveling Salesman Pairs
         for ifile = [1 fidtxtout]
            fprintf(ifile,'%3d %3d %3d %3d %5d %5d %12.4f %12.4f %12.4f %12.4f\n',i,j,k,k+1 ...
               ,iuniqorbs(i0),iuniqorbs(i1),tepochs(i0),  tepochs(i1)...
               ,tepochs(i1)-tepochs(i0)...
               ,bpest(i1)-  bpest(i0));
         end
         % connect TSP with dotted line
         if plotts < 3
            plot([tepochs(i0) tepochs(i1)],[bpest(i0) bpest(i1)]...
               ,mylins{1+mod(j,length(mylins))},'Linewidth',2);
         end
      end

   end
end


for i=1:np
   for j = 1:nf
      if sum(ismember(species(j,:),id0(i))) == 1 && sum(ismember(species(j,:),id1(i))) == 1
         % draw symbol
%          plot([tepochs(id0(i)) tepochs(id1(i))],[bpest(id0(i)) bpest(id1(i))],mysym0{1+mod(j,length(mysym0))},'Linewidth',2,'MarkerFaceColor','k'); hold on;
        plot([tepochs(id0(i)) tepochs(id1(i))],[bpest(id0(i)) bpest(id1(i))],mysym0{sat_id(i)},'Linewidth',2,'MarkerSize', 10); hold on;
          
         if plotts == 0  || plotts == 2 || plotts == 3 % draw dashed line for possible pairs
%               plot([tepochs(id0(i)) tepochs(id1(i))],[bpest(id0(i)) bpest(id1(i))],mydash{1+mod(j,length(mydash))},'Linewidth',2,'MarkerFaceColor','k'); hold on;
              plot([tepochs(id0(i)) tepochs(id1(i))],[bpest(id0(i)) bpest(id1(i))],mydash{sat_id(i)},'Linewidth',2,'MarkerSize', 10); hold on;

         end
      end
   end
end



plot([min(tepochs) min(tepochs)],[min(bpest)-0.1*(max(bpest)-min(bpest)) max(bpest)+0.1*(max(bpest)-min(bpest))],'w.'); % draw a white dot to stretch scales


v = axis;

if lab_on == 1
% Plotting Epoch Calendar Dates
kt = floor(min(tepochs)):1:ceil(max(tepochs)); %getting time span

mepi = 0; %set max number of epochs per interval = 0 
for k = 1:length(kt)-1
    Icount = find(tepochs > kt(k) & tepochs <= kt(k+1)); % set vector of years 
    if length(Icount) > mepi
        mepi = length(Icount); %finding maximum number of epochs per interval
    end
end

maxepi = ceil(mepi/3); 
deltay = ((max(bpest)-min(bpest))/5); % find distance along y axis for writing
maxy = min(bpest)+deltay*7; % counts 5 deltay for length of longest bpest plus 2 for top border
miny = min(bpest)-deltay*maxepi; % sets bottom border for text

month_name = {'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'}; % string of month names
[yr mn day] = dyear2date(tepochs);

del = deltay/2; %amount of space between each label, dependent on tick marks 

lspace = min(bpest)-miny;%ceil(min(bpest)-ceil(mepi*pix)); %finding amount of space needed for labels
l_bot = miny;%lspace+del; %starting spot for label
for k = 1:length(kt)-1
    I = find(yr >= kt(k) & yr < kt(k+1));
    xp = kt(k);%(kt(k)+kt(k+1))/2; %average the x interval for text starting point
    tp = yr(I); 
    dy = length(tp);
    by = l_bot +del*(dy-1); %9.5-.75*(dy-1); -14.5
    yy = fliplr(l_bot:del:by); %fliplr(by:.75:9.5);
    for xx = 1:dy
       % text(xp,yy(xx),sprintf('%1s-%1d',month_name{mn(I(xx))}, day(I(xx)))); % print calendar month and day for each epoch in interval
       text(xp,yy(xx),sprintf('%3s-%02.2d',month_name{mn(I(xx))}, day(I(xx)))); % print calendar month and day for each epoch in interval
        
    end
end
axis([floor(min(tepochs))-1 ceil(max(tepochs))+1 miny-deltay maxy])
end

% Set axis and label properties
year_counts = unique(floor(tepochs));
year_lab = floor(min(tepochs)):2:ceil(max(tepochs)); % set x axis labels
year_labs = {};
for i = 1:numel(year_lab)
    j = i+1*(i-1);
    year_labs{j} = num2str(year_lab(i));
    year_labs{j+1} = ' ';
end
set(gca,'XTick',[floor(min(tepochs)):1:ceil(max(tepochs))]) 
set(gca,'XTickLabel',year_labs)

% Set figure properties
h2=title ('Minimum Spanning Forest (MSF) pairs'); set(h2,'FontName','Courier','Fontsize',14,'FontWeight','bold');
h2=xlabel(xlab);       set(h2,'FontName','Courier','Fontsize',14,'FontWeight','bold');
% h2=ylabel(ylab,'interpreter','latex');         set(h2,'FontName','Courier','Fontsize',14,'FontWeight','bold');
h2=ylabel('volumetric strain rates [1/yr]');         set(h2,'FontName','Courier','Fontsize',14,'FontWeight','bold');
set(gca,'FontName','Courier','Fontsize',14,'FontWeight','bold');

if nargout == 0
   hold off
end




for ifile = [1 fidtxtout]
   fprintf(ifile,'%s ended at %s\n',mfilename,datestr(now,31));
end
fclose(fidtxtout);

 set(FigHandle, 'Position', [100, 100, 995, 895]);%1049
%save('OkmokBplot')
return;


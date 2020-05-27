function plot_trees3(tepochs, yvals, Qiev, itrees, xlab, ylab, titlestr)
%function plot_trees3(tepochs, yvals, Qiev, trees,xlab, ylab)
%
% plot graph 
% for example pseudo-absolute Bperp as a function of time 
% inputs:
%     tepochs == decimal years corresponding to vertices (epochs)
%     scores  == weighting of edges
%     itrees   == list of indices to vertices in trees, one row per tree
%     Qiev    == Edge-vertex incidence matrix
% outputs: 
%     Hfig    == graphics handle of figure
%
% Automatically handle 3 cases
% 1  X = t      Y = Bperp
% 2  X = t      Y = Ddop
% 3  X = Bperp  Y = Ddop
%
%
% ktours = plot_trees(tepochs, scores, DD, trees, iuniqorbs, uniqdates, plotts)
%
%
% ktours = plot_trees(tepochs, scores, DD, trees)
%  
% ktours = plot_trees(tepochs, scores, DD, trees, iuniqorbs, uniqdates)
%
% ktours = plot_trees(tepochs, scores, DD, trees, iuniqorbs, uniqdates, plotts)
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
% 2020-MAR-26 updated comments
% 2020-APR-24 do not return graphics handle

fidtxtout = fopen(sprintf('%s_out.txt',mfilename),'a+t');
% for ifile = [1 fidtxtout]
%    fprintf(ifile,'%s begins at %s\n',mfilename,datestr(now,31));
% end


narginchk(4,8);
if nargin <= 6
   plotts = 0;
end
nargoutchk(0,0);

if exist('xlab','var') == 0
    xlab = 'year';
end
if exist('ylab','var') == 0
    xlab = 'relative score';
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

mysyms = {'gx-' 'ro-' 'b*-' 'ks-' 'md-' 'cv-'};
mysols = {'g-'  'r-'  'b-'  'k-'  'm-'  'c-'};
mylins = {'g:'  'r:'  'b:'  'k:'  'm:'  'c:'};
mydash = {'g--' 'r--' 'b--' 'k--' 'm--' 'c--'};
mysym0 = {'gx'  'ro'  'b*'  'ks'  'md'  'cv'};

% graphics handle to return
%h=figure;hold on; 

% number of trees
[ntrees,ndum] = size(itrees);

% number of pairs and number of epochs
[np,me] = size(Qiev); 

%% Kurt 20150901 Causes problems because iuniqorbs not defined
% %plot origin to make legend come out right
% if ntrees < 10
%     for j=1:ntrees
%         plot(min(tepochs),0,mysyms{1+mod(j,length(mysyms))});
%         tree = trees(j,:);
%         k=isfinite(tree);
%         tree=tree(k);
%         me = length(find(k == 1));
%         if nargin >= 6
%             tree_name{j} = strcat(sprintf('trees %s orbits:',char(j+64)),sprintf(' %7d',iuniqorbs(tree(1:me))));
%         else
%             tree_name{j} = strcat(sprintf('trees %s ID:',char(j+64)),sprintf('%3d',char(j+64),tree(1:me)));
%         end
%     end
%     ktours = zeros(size(trees));
%     
%     legend(tree_name,'Location','NorthOutside');
%     % over plot origin with white
%     plot(min(tepochs),0,'sw');
%     plot(min(tepochs),0,'ow');
%     plot(min(tepochs),0,'xw');
%     plot(min(tepochs),0,'*w');
%     plot(min(tepochs),0,'dw');
%     plot(min(tepochs),0,'vw');
%     
% end


% %% draw end points of available pairs
% id0 = zeros(np,1); 
% id1 = zeros(np,1); 
% for i=1:np      
%    ddcol = Qiev(i,:); 
%    j=find(abs(ddcol)>0);
%    id0(i) = min(j);  % index to first  (master) epoch
%    id1(i) = max(j);  % index to second (slave) epoch
% end

iedges = incidence_to_edges(Qiev);
id0 = iedges(:,1); % index to starting vertex in graph, corresponding to (master) epoch
id1 = iedges(:,2); % index to finishing  vertex in graph corresponding to (slave) epoch
figure;
hold on;
for i=1:np
   for j = 1:ntrees
      if sum(ismember(itrees(j,:),id0(i))) == 1 && sum(ismember(itrees(j,:),id1(i))) == 1  
         % draw symbol at vertices
         plot([tepochs(id0(i)) tepochs(id1(i))],[yvals(id0(i)) yvals(id1(i))],mysym0{1+mod(j,length(mysym0))},'Linewidth',2,'MarkerFaceColor','k'); hold on;
         
         % draw dashed line on edges (pairs)  
         plot([tepochs(id0(i)) tepochs(id1(i))],[yvals(id0(i)) yvals(id1(i))],mydash{1+mod(j,length(mydash))},'Linewidth',2,'MarkerFaceColor','k'); hold on;
      end
   end
end




%plot([min(tepochs) min(tepochs)],[min(scores)-0.1*(max(scores)-min(scores)) max(scores)+0.1*(max(scores)-min(scores))],'w.'); % draw a white dot to stretch scales

% %% Kurt 2015 causes problems
% % make a legend
% if np < 10
%    legend(tree_name,'Location','NorthOutside');
% end

% make title string
titl = strcat(titlestr,': ',sprintf('np = %d me = %d ntrees = %d\n',np,me,ntrees));

title (titl,'FontName','Helvetica','Fontsize',14,'FontWeight','bold');
xlabel(xlab,'FontName','Helvetica','Fontsize',14,'FontWeight','bold');
ylabel(ylab,'FontName','Helvetica','Fontsize',14,'FontWeight','bold');
set(gca,'FontName','Helvetica','Fontsize',14,'FontWeight','bold');

% for ifile = [1 fidtxtout]
%    fprintf(ifile,'%s ended at %s\n',mfilename,datestr(now,31));
% end
fclose(fidtxtout);

% % return graphics handle
% Hfig = gcf;
return;
end


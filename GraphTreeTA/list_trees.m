function treeList = list_trees(itrees,vertexIds)
% % given a list of trees indexing vertices, make a list of vertices
% %
% 
%     %Trees=cell(ntrees,1);
%     itrees = nan(ntrees,nvertices);
%     ivertices = 1:nvertices;
%     for j=1:ntrees
%         inonnull = find(abs(N(:,j)) > 0);
%         % get row vector of indices of vertices
%         ivertices1 = rowvec(ivertices(inonnull));
%         % store indices of epochs in jth component
%         %Trees{j} = ivertices1; 
%         n1 = numel(ivertices1);
%         itrees(j,1:n1) = ivertices1;
%         for k=n1+1:nvertices
%            itrees(j,k)=nan;
%         end
%     end
% end
% 
% %Trees = nan;
% 
% % % make one tree per row
% % Trees = Trees';
% % 
% % % print out
% % if verbose == 1   
% %     for i=1:ntrees
% %         fprintf(1,'Tree%03d contains vertices: ',i);
% %         fprintf(1,'%d ',Trees{i});
% %         fprintf(1,'\n');
% %     end
% % end
% % 
% % 
% % % sort trees
% % %Trees = sortrows(Trees);
% end


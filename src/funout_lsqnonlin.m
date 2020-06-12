function [state,options,optchanged] = funout_lsqnonlin(options,state,flag)
% persistent h1 history r
% optchanged = false;
switch flag
    case 'init'
        %        h1 = figure;
        %         ax = gca;
        %         ax.XLim = [0 21];
        %         ax.YLim = [0 21];
        %         l1 = min(state.Population(:,1));
        %         m1 = max(state.Population(:,1));
        %         l2 = min(state.Population(:,2));
        %         m2 = max(state.Population(:,2));
        %         r = rectangle(ax,'Position',[l1 l2 m1-l1 m2-l2]);
        %         history(:,:,1) = state.Population;
        %         assignin('base','gapopulationhistory',history);
        fprintf(1,'%s init\n',mfilename);
    case 'iter'
        %         % Update the history every 10 generations.
        %         if rem(state.Generation,10) == 0
        %             ss = size(history,3);
        %             history(:,:,ss+1) = state.Population;
        %             assignin('base','gapopulationhistory',history);
        %         end
        %         % Find the best objective function, and stop if it is low.
        %         ibest = state.Best(end);
        %         ibest = find(state.Score == ibest,1,'last');
        %         bestx = state.Population(ibest,:);
        %         bestf = gaintobj(bestx);
        %         if bestf <= 0.1
        %             state.StopFlag = 'y';
        %             disp('Got below 0.1')
        %         end
        %         % Update the plot.
        %         figure(h1)
        %         l1 = min(state.Population(:,1));
        %         m1 = max(state.Population(:,1));
        %         l2 = min(state.Population(:,2));
        %         m2 = max(state.Population(:,2));
        %         r.Position = [l1 l2 m1-l1 m2-l2];
        %         pause(0.1)
        %         % Update the fraction of mutation and crossover after 25 generations.
        %         if state.Generation == 25
        %             options.CrossoverFraction = 0.8;
        %             optchanged = true;
        %         end
        fprintf(1,'%s iter\n',mfilename);
        
    case 'done'
        %         % Include the final population in the history.
        %         ss = size(history,3);
        %         history(:,:,ss+1) = state.Population;
        %         assignin('base','gapopulationhistory',history);
        fprintf(1,'%s done\n',mfilename);
        
end
return
end
function ok=fexist(fname)
% % check to see if a file named fname exists

% home-made version
ok=fopen(fname,'r');
if ok > 0
   fclose(ok);
   ok = 1;
else
   ok = 0;
end

% Use built-in function instead
% % This looks anywhere in the Matlab command search path
% if exist(fname,'file')  == 2
%     ok = 1;
% else
%     ok = 0;
% end
return

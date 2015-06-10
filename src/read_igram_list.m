function [mdate, imast, sdate, islav, hamb, ddays, t0, t1] = read_igram_list(igram_list_file,key)
% read a list of interferometric pairs, 1 per line, with 'key' on the line
% [mdate, imast, sdate, islav, hamb, ddays, t0, t1] = read_igram_list(igram_list_file,key)
% for use with DIAPASON DTOOLS
fprintf(1,'%s begins ...\n',mfilename);

fid = fopen(igram_list_file,'r');
if fid <= 0 
    error 'Cannot open list of interferograms'
    return
end
i = 0;
while 1 % for each pair
    line1 = fgetl(fid);
    if ~isstr(line1), break, end
    if findstr(line1,key)
        i = i+1;
        ftmp=fopen('tmp.txt','w'); % write a scratch file
        fprintf(ftmp,'%s',line1);
        fclose(ftmp);
        %     [myr(i) uobso{i} mdy(i) syr(i) smo{i} sdy(i) imast(i) islav(i) hamb(i) ddays(i) t0(i) t1(i)] = textread ('adj.tmp','%d %3s %d %d %3s %d %d %d %f %f %f %f%*[^\n]');
        [myr mmo mdy syr smo sdy imast(i) islav(i) hamb(i) ddays(i) t0(i) t1(i)] = textread ('tmp.txt','%s %s %s %s %s %s %d %d %f %f %f %f%*[^\n]');
        %     mstr = sprintf('%s-%s-%s',myr{1},mmo{1},mdy{1});
        %     sstr = sprintf('%s-%s-%s',syr{1},smo{1},sdy{1});
        %     2009-DEC-02 pad with zeros
        if numel(str2num(myr{1}))>0 & numel(str2num(mmo{1}))>0 & numel(str2num(mdy{1}))>0
            mstr = sprintf('%04d-%02d-%02d',str2num(myr{1}),str2num(mmo{1}),str2num(mdy{1}));
        else
            mstr = sprintf('%s-%s-%s',myr{1},mmo{1},mdy{1});
        end
        if numel(str2num(syr{1}))>0 & numel(str2num(smo{1}))>0 & numel(str2num(sdy{1}))>0
            sstr = sprintf('%04d-%02d-%02d',str2num(syr{1}),str2num(smo{1}),str2num(sdy{1}));
        else
            sstr = sprintf('%s-%s-%s',syr{1},smo{1},sdy{1});
        end
        
        mdate{i} = mstr;
        sdate{i} = sstr;
    end
end

if (i == 0) 
    fprintf(1,'ERROR No lines in %s matched key %s\n',igram_list_file,key);
    error 'No records matched key' 
    return
end
% convert to column vectors
%mdate = mdate';
%sdate = sdate';
imast = imast';
islav = islav';
hamb = hamb';
ddays = ddays';
t0 = t0';
t1 = t1';

fprintf(1,'Read %4d pairs containing key %s\n',i,key);

fclose(fid);
return


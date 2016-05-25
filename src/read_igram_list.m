function [mdate, imast, sdate, islav, hamb, ddays, t1, t2] = read_igram_list(igram_list_file,key)
% read a list of interferometric pairs, 1 per line, with 'key' on the line
% [mdate, imast, sdate, islav, hamb, ddays, t0, t1] = read_igram_list(igram_list_file,key)
% for use with DIAPASON DTOOLS
% 20160524 use datetime structures for t1 and t2
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
        %[myr(i) mmostr{i} mdy(i) syr(i) smostr{i} sdy(i) imast(i) islav(i) hamb(i) ddays(i) y1(i) y2(i)] = textread('tmp.txt','%d %s %d %d %s %d %d %d %f %f %f %f%*[^\n]')
        [myr mmostr mdy syr smostr sdy imast(i) islav(i) hamb(i) ddays0 y1 y2] = textread ('tmp.txt','%d %s %d %d %s %d %d %d %f %f %f %f%*[^\n]')
        %     mstr = sprintf('%s-%s-%s',myr{1},mmo{1},mdy{1});
        %     sstr = sprintf('%s-%s-%s',syr{1},smo{1},sdy{1});
        %     2009-DEC-02 pad with zeros

        mmo = monthstr2monthnum(char(mmostr{1}));
        tm = datetime(myr,mmo,mdy,'TimeZone','UTC');
        tm.Format = 'yyyy-MM-dd'  ; % 24 hour clock
        mstr = sprintf('%s',tm);
        
        smo = monthstr2monthnum(char(smostr{1}));
        ts = datetime(syr,smo,sdy,'TimeZone','UTC');
        ts.Format = 'yyyy-MM-dd'  ; % 24 hour clock
        sstr = sprintf('%s',ts);
         
        mdate{i} = mstr;
        sdate{i} = sstr;
        
        t1(i) = tm;
        t2(i) = ts;
        ddays(i) = days(ts-tm);
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
% t1 = t1';
% t2 = t2';

fprintf(1,'Read %4d pairs containing key %s\n',i,key);

fclose(fid);
return


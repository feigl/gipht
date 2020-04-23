function Tgps = get_gps(site) 
%% get GPS data for station named 'site'


% http://geodesy.unr.edu/gps_timeseries/tenv3/IGS14/GARL.tenv3

%! curl http://geodesy.unr.edu/gps_timeseries/tenv3/IGS14/GARL.tenv3 -o /Users/feigl/BoxSync/WHOLESCALE/GPS/GARL.tenv3

urlFolder='http://geodesy.unr.edu/gps_timeseries/tenv3/IGS14';
%gpsFileName1 = 'GARL.tenv3';
gpsFileName1 = strtrim(sprintf('%s.tenv3',site));

%gpsFileName2 = websave('GARL.tenv3','http://geodesy.unr.edu/gps_timeseries/tenv3/IGS14/GARL.tenv3');
gpsFileName2 = websave(gpsFileName1,strcat(urlFolder,'/',gpsFileName1));

nHeaderLines = 1;
Tgps = readFlatFileAsTable(gpsFileName1,nHeaderLines);
[ndata, ncols] = size(Tgps);

fprintf(1,'\nFirst 10 lines:\n');
Tgps(1:10,{'yyyy_yyyy','u0_m_','up_m_'})
fprintf(1,'Last 10 lines:\n');
Tgps(end-10:end,{'yyyy_yyyy','u0_m_','up_m_'})
return
end






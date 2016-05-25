function ierr=read4(fname,ncols,dirname)
%function ierr=read4(fname,ncols,dirname)
% read and translate a file from OKMOK
% Kurt 2009-JAN-16
% modified into a function by Summer 2009-APR-01
% 2009-05-08 Kurt
% 2011-06-12 Update for new,new 20-meter interferograms from Zhong Lu
%
% Try two ways, one (or both) of which will solve the problem.
% 1) Treat each 4-byte part (i.e., the real or imaginary part of the fn{1} image, or the am{1} image) as a float-type data. If it still does not solve the problem,
% 2) Please swap the bytes when you read in these images.
% hengill.geology.wisc.edu% dir ../../old_okmok/238274_238775*
% -rw-rw-rw-  1 summer ohlendor  385 Mar 11  2008 ../../old_okmok/238274_238775.base
% -rw-rw-rw-  1 summer ohlendor 4.4M Mar 11  2008 ../../old_okmok/238274_238775.pwr1.utm1
% -rw-rw-rw-  1 summer ohlendor 8.7M Mar 11  2008
% ../../old_okmok/238274_238775.sub_int.utm1
% Example:
%   ierr=read4('227524_228526.sub_int.utm1',2200,'../E2_T115');

%         98765432109876543210
ierr = 1;

%ip=findstr(fname,'.pwr1')
ip=findstr(fname,'.utm');
ip=findstr(fname,'.sub_int.utm1');
e15=str2num(fname(ip-11:ip-7)); %changed from e1,e2 that were for long orbit #s with sat code
e25=str2num(fname(ip-5 :ip-1));

% remove leading digit coding satellite NO SAT CODE FOR ENV
%e15 = mod(e1,100000);
%e25 = mod(e2,100000);
% make long path name for file containing phase and coherence
%pname = sprintf('../old_okmok/%d_%d.sub_int.utm',e1,e2)
% make long path name for file containing amplitude
%aname = sprintf('../old_okmok/%d_%d.pwr1.utm1',e1,e2)

if (e15<10000)
	if (e25<10000)
		pname = sprintf('%s/0%d_0%d.sub_int.utm1',dirname,e15,e25);
		%aname = sprintf('%s/0%d_0%d.pwr1.utm1',dirname,e15,e25); %don't
		%have these for ENV T115 and 344 from ZL- use master's amp
        aname = sprintf('%s/0%d.pwr.utm1',dirname,e15);  
	else
		pname = sprintf('%s/0%d_%d.sub_int.utm1',dirname,e15,e25);
        %aname = sprintf('%s/0%d_%d.pwr1.utm1',dirname,e15,e25);
        aname = sprintf('%s/0%d.pwr.utm1',dirname,e15);
	end
else
	if (e25<10000)
		pname = sprintf('%s/%d_0%d.sub_int.utm1',dirname,e15,e25);
        %aname = sprintf('%s/%d_0%d.pwr1.utm1',dirname,e15,e25);
        aname = sprintf('%s/%d.pwr.utm1',dirname,e15);
	else
		pname = sprintf('%s/%d_%d.sub_int.utm1',dirname,e15,e25);
        %aname = sprintf('%s/%d_%d.pwr1.utm1',dirname,e15,e25);
        aname = sprintf('%s/%d.pwr.utm1',dirname,e15);
	end
end

z4=read_swap_cr4(pname,ncols);
% changed 04/14/09 to reverse phase progression
%pha = angle(z4) / pi / 2; % from -0.5 to +0.5 cycles
pha = -1* angle(z4) / pi / 2; % from -0.5 to +0.5 cycles such that increasing phase is increasing range
coh = abs(z4)/(max(max(abs(z4))-min(min(z4)))); % from 0 to 1
ibad=find(abs(coh) ==0);
fprintf(1,'Number of bad pixels in phase image = %ld\n',numel(ibad));
pha(ibad) = NaN;
coh(ibad) = NaN;

% Read the Amplitude file
if fexist(aname) == 0
    amp = ones(size(pha));
else
    amp=read_swap_r4(aname,ncols);
end

% Mask pixels that are bad in phase or coherence
%amp(ibad) = NaN;

% mask pixels with zero amplitude
% ibad=find(amp == 0);
% amp(ibad) = NaN;

% log stretch
amp=100*log(amp);
ibad=find(isfinite(amp) == 0);
amp(ibad)=NaN;
% truncat bottom 5 percent
ampmin=quantile(reshape(amp,numel(amp),1),0.05)
imin=find(amp < ampmin & isfinite(amp)==1);
amp(imin)=ampmin;

% truncate top 5 percent
ampmax=quantile(reshape(amp,numel(amp),1),0.95)
imax=find(amp > ampmax & isfinite(amp) ==1);
amp(imax)=ampmax;

% rescale to range from 0 to 1
amp = (amp - ampmin) / (ampmax - ampmin);
ampmax=max(max(amp))
ampmin=min(min(amp))


figure
colormap('gray');
%colormap('jet');
subplot(3,1,1);
imagesc(pha);hold on; colorbar;title('PSP phase in cycles');axis image;
subplot(3,1,2);
imagesc(coh);hold on; colorbar;title('coherence'); axis image
subplot(3,1,3);
imagesc(amp);hold on; colorbar;title('amplitude'); axis image
printjpg('3panelraw.jpg')

figure
%    subplot(2,2,1);hist(reshape(real(z4),numel(z4),1),40);title('real part');
%    subplot(2,2,2);hist(reshape(imag(z4),numel(z4),1),40);title('imag part');
subplot(2,2,1);hist2(reshape(abs(z4),numel(z4) ,1),40);title('absolute value');
subplot(2,2,2);hist2(reshape(amp,    numel(amp),1),40);title('amplitude');
subplot(2,2,3);hist2(reshape(pha,    numel(pha),1),40);title('phase in cycles');
subplot(2,2,4);hist2(reshape(coh,    numel(coh),1),40);title('coherence');
printjpg('hists.jpg')

% phafile=sprintf('../../old_okmok/pha_%d_%d_ort.pha',e1,e2)
% pspfile=sprintf('../../old_okmok/smp_%d_%d_ort.pha',e1,e2)
% cohfile=sprintf('../../old_okmok/coh_%d_%d_ort.byt',e1,e2)
% ampfile=sprintf('../../old_okmok/amp_%d_%d_ort.byt',e1,e2)
phafile=sprintf('%s/pha_%d_%d_ort.pha',dirname,e15,e25)
pspfile=sprintf('%s/psp_%d_%d_ort.pha',dirname,e15,e25)
cohfile=sprintf('%s/coh_%d_%d_ort.byt',dirname,e15,e25)
ampfile=sprintf('%s/amp_%d_%d_ort.byt',dirname,e15,e25)

write_pha(phafile,floor(pha*256));
write_byt(cohfile,floor(coh*256));
write_byt(ampfile,floor(amp*256));

% Power-spectrum filtering using unit amplitude
%pscmd = sprintf('/usr/local/diapason4/utils/ps_filt2 %s ones.byt %s %d 0.7\n',phafile,pspfile,ncols)
% Power spectrum filtering using observed amplitude
pscmd = sprintf('/usr/local/diapason4/utils/ps_filt2 %s %s %s %d 0.7\n',phafile,ampfile,pspfile,ncols)
[s, w] = unix(pscmd);

% check what we have written
amp2=read_byt(ampfile,ncols);
smp2=read_pha(pspfile,ncols);
coh2=read_byt(cohfile,ncols);

figure
subplot(3,1,1); colormap('jet');
imagesc(smp2);hold on; colorbar;title('SMP phase [-128,+127]');axis image;
subplot(3,1,2);  colormap('copper');
imagesc(coh2);hold on; colorbar;title('coherence [0,255]'); axis image;
subplot(3,1,3);    colormap('gray');
imagesc(amp2);hold on; colorbar;title('amplitude [0,255]'); axis image;
printjpg('3panelsmp.jpg');

ierr = 0;
return






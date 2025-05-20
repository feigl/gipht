function [V2,tstring,iref,jref,kref,kok,V0] = extract_and_reference_velocity(AV,VGRD,SGRD,XGRD,YGRD,IGRD,UTMranges,iReferencePoint,Eref,Nref)
%function V2 = extract_and_reference_velocity(DV,AV,UTMranges)
% 2024/09/28 Kurt Feigl
% 2025/04/10 added UTM E,N for reference pixel

kilo=1000;
% use logicals to find pixels within bounds
iok=and((XGRD >= UTMranges(1)),(XGRD <= UTMranges(2)));%size(iok)
jok=and((YGRD >= UTMranges(3)),(YGRD <= UTMranges(4)));%size(jok)
kok=and(iok,jok);

% subset on finite values
kok=and(kok,isfinite(VGRD));

% subset pixels with SNR > 2
if isfinite(SGRD) == true
   whos *GRD
   kok=and(kok,isfinite(VGRD));
   kok=and(kok,(abs(VGRD ./ SGRD) > 2.));
else
    warning('SGRD is not defined.');
end

%iReferencePoint=3
switch iReferencePoint
    case 1
        if isfinite(VGRD(AV.REF_X,AV.REF_Y))
            tstring = sprintf('w.r.t. most coherent pixel');
            V0=VGRD(AV.REF_X,AV.REF_Y);
            kref=sub2ind(size(VGRD),AV.REF_X,AV.REF_Y);
        else
            warning("velocity at most coherent pixel is NaN")
            V0=zeros(size(VGRD));
            tstring = sprintf('from MintPy');
            [~,kref]=find(isfinite(VGRD), 1 );
        end
        [iref,jref]=ind2sub(size(VGRD),kref);
    case 2
        [~,kref]=min(SGRD,[],'all');
        [iref,jref]=ind2sub(size(SGRD),kref);
        V0=VGRD(kref);
        tstring = sprintf('w.r.t. pixel with minimum std dev');
    case 3
        % use logicals for SW corner
        iref=and((XGRD - UTMranges(1) <= 2000), (XGRD - UTMranges(1) >0));
        jref=and((YGRD - UTMranges(3) <= 2000), (YGRD - UTMranges(3) >0));
        kref=and(jref,iref);
        V0=median(VGRD(kref),'all','omitnan');
        tstring = sprintf('w.r.t. median SW corner' );
    case 4
        % use logicals for area around reference
        iref=and(iok,(abs(XGRD - Eref) <= 1000));
        jref=and(jok,(abs(YGRD - Nref) <= 1000));
        kref=and(jref,iref);
        V0=median(VGRD(kref),'all','omitnan');
        tstring = sprintf('w.r.t. median square at (E,N) =(%.3f %.3f) [km]',Eref/kilo,Nref/kilo);
   case 5
        % use logicals for NW corner
        iref=and((XGRD - UTMranges(1) <= 2000), (XGRD - UTMranges(1) >     0));
        jref=and((YGRD - UTMranges(4) <=   0),  (YGRD - UTMranges(4) >= -2000));
        kref=and(jref,iref);
        V0=median(VGRD(kref),'all','omitnan');
        tstring = sprintf('w.r.t. SW corner' );
   case 6
        % use logicals for NE corner
        iref=and((XGRD - UTMranges(2) <=   0),  (XGRD - UTMranges(2) >= -2000));
        jref=and((YGRD - UTMranges(4) <=   0),  (YGRD - UTMranges(4) >= -2000));
        kref=and(jref,iref);
        V0=median(VGRD(kref),'all','omitnan');
        tstring = sprintf('w.r.t. NE corner' );
otherwise
        error('unknown iReferencePoint %d',iReferencePoint)
end


V2=double((VGRD - V0) ./ cosd(IGRD));

kok=and(kok,isfinite(kok));

nok=numel((kok))
if nok < 100
    warning('Very few points %d',nok);
end

return
end


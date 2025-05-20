function [V2,tstring,iref,jref,kref,kok,V0] = extract_and_reference_velocity(DV,AV,UTMranges,iReferencePoint,Eref,Nref)
%function V2 = extract_and_reference_velocity(DV,AV,UTMranges)
% 2024/09/28 Kurt Feigl
% 2025/04/10 added UTM E,N for reference pixel

kilo=1000;
% use logicals to find pixels within bounds
iok=and((DV.XGRD >= UTMranges(1)),(DV.XGRD <= UTMranges(2)));%size(iok)
jok=and((DV.YGRD >= UTMranges(3)),(DV.YGRD <= UTMranges(4)));%size(jok)
kok=and(iok,jok);
% subset on finite values
kok=and(kok,isfinite(DV.velocity));

% subset pixels with SNR > 2
kok=and(kok,(abs(DV.velocity ./ DV.velocityStd) > 2.));

%iReferencePoint=3
switch iReferencePoint
    case 1
        if isfinite(DV.velocity(AV.REF_X,AV.REF_Y))
            tstring = sprintf('Mean LOS velocity w.r.t. most coherent pixel [mm/year]');
            V0=DV.velocity(AV.REF_X,AV.REF_Y);
            kref=sub2ind(size(DV.velocity),AV.REF_X,AV.REF_Y);
        else
            warning("velocity at most coherent pixel is NaN")
            V0=zeros(size(DV.velocity));
            tstring = sprintf('Mean LOS velocity from MintPy');
            [~,kref]=find(isfinite(DV.velocity), 1 );
        end
        [iref,jref]=ind2sub(size(DV.velocity),kref);
    case 2
        [~,kref]=min(DV.velocityStd,[],'all');
        [iref,jref]=ind2sub(size(DV.velocityStd),kref);
        V0=DV.velocity(kref);
        tstring = sprintf('Mean LOS velocity w.r.t. pixel with lowest std [mm/year]');
    case 3
        % use logicals for SW corner
        iref=and((DV.XGRD - UTMranges(1) <= 400), (DV.XGRD - UTMranges(1) >0));
        jref=and((DV.YGRD - UTMranges(3) <= 400), (DV.YGRD - UTMranges(3) >0));
        kref=and(jref,iref);
        V0=mean(DV.velocity(kref),'all','omitnan');
        tstring = sprintf('Mean LOS velocity w.r.t. SW corner [mm/year]' );
    case 4
        % use logicals for area around reference
        iref=and(iok,(abs(DV.XGRD - Eref) <= 200));
        jref=and(jok,(abs(DV.YGRD - Nref) <= 200));
        kref=and(jref,iref);
        V0=mean(DV.velocity(kref),'all','omitnan'));
        tstring = sprintf('Mean LOS velocity [mm/year] w.r.t. square at (E,N) =(%.3f %.3f) [km',Eref/kilo,Nref/kilo);
    otherwise
        error('unknown iReferencePoint %d',iReferencePoint)
end

V2=double(DV.velocity - V0);

tstring = strcat(tstring,sprintf(' %8d to %8d',AV.START_DATE,AV.END_DATE));
return
end


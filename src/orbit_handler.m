function [unitv, orbvm, orbvs] = orbit_handler(i, i1, i2, orbfile, orbits_loaded, orbopt, unitv0)
%function [unitv] = orbit_handler(i, i1, i2, orbfile, orbits_loaded, orbopt, unitv0)
% only tested for orbfile = ''
    if orbits_loaded == 0
        % orbital information varies across scene and thus across pairs
        unitv(1:3,i1:i2) = nan;
        orbvm(1:6,i1:i2) = 0.0; % partial derivative with respect to master orbital parameters
        orbvs(1:6,i1:i2) = 0.0; % partial derivative with respect to slave  orbital parameters
        if numel(orbfile) > 0
            % Calculate orbital quantities at center for master
            [xs, ys, zs, xdot, ydot, zdot, mjd, secs, orbnum] = readorb(ofnames1{i});
            if min(orbnum) ~= abs(imast(i))
                warning(sprintf('Orbit numbers differ %d %d %d\n',min(orbnum),abs(imast(i)),min(orbnum)-abs(imast(i))));
            end
            [UC, DC, NC, RC, HC, AC, VC, near_mjdC, near_secC] = orbvectors2(loncenter,latcenter,0.0 ...
                ,xs, ys, zs, xdot, ydot, zdot, mjd, secs, orbnum);
            SC = DC/norm(DC); % unit look vector target w.r.t. satellite
            
            for ii=i1:i2
                if ii == i1
                    fprintf(1,'\nCalculating look vectors for MASTER of pair %d for pixels %d to %d from orbit file %s\nPixels numbered:\n',i,i1,i2,ofnames1{i});
                    tstart=tic;
                end
                % calculate unit vector that varies across scene
                if mod(ii,10) == 0; fprintf(1,'%4d ',ii); end
                if mod(ii,100) == 0; fprintf(1,'\n'); end
                %  calculate components of orbit vector
                [U1, D1, N1, R1, H1, A1, V1, near_mjd1, near_sec1] = orbvectors2(alon(ii),alat(ii),xyzm(3,ii)...
                    ,xs, ys, zs, xdot, ydot, zdot, mjd, secs, orbnum);
                if abs(near_mjd1 - near_mjdC) > 0
                    fprintf(1,'Longitudes %12.6f %12.6f %12.6f\n',loncenter,alon(ii),loncenter-alon(ii));
                    fprintf(1,'Latitudes  %12.6f %12.6f %12.6f\n',latcenter,alat(ii),latcenter-alat(ii));
                    dt = 86400.0*(near_mjd1 - near_mjdC) + (near_sec1 - near_secC); % time difference in seconds
                    warning(sprintf('different MJD %d %d %d\n',near_mjd1,near_mjdC,near_mjd1 - near_mjdC));
                else
                    dt = near_sec1 - near_secC; % time difference in seconds
                end
                
                unitv(1,ii) = U1(1);    % eastward  component of dimensionless unit vector
                unitv(2,ii) = U1(2);    % northward component of dimensionless unit vector
                unitv(3,ii) = U1(3);    % upward    component of dimensionless unit vector
                
                if orbopt == 1
                    SC = D1/norm(D1); % unit look vector target w.r.t. satellite
                    
                    if ismember(pselect,[7,9]) == 1
                        %  calculate components of baseline
                        [U3, D3, N3, R3, H3, A3, V3, near_mjd3, near_sec3] = orbvectors2(alon(ii)+dlon,alat(ii),xyzm(3,ii)...
                            ,xs, ys, zs, xdot, ydot, zdot, mjd, secs, orbnum);
                        S3 = D3/norm(D3); % unit look vector target w.r.t. satellite
                        orbvm(1,ii) = dot(H3,SC) - dot(H1,SC); % partial derivative of range w.r.t. horizontal  component of orbit adjustment vector
                        orbvm(2,ii) = dot(A3,SC) - dot(A1,SC); % partial derivative of range w.r.t. along-track component of orbit adjustment vector
                        orbvm(3,ii) = dot(V3,SC)-  dot(V1,SC); % partial derivative of range w.r.t. vertical    component of orbit adjustment vector
                        orbvm(4,ii) = dt * (dot(H3,SC)-dot(H1,SC)); % partial derivative of range w.r.t. horizontal  component of orbit adjustment velocity
                        orbvm(5,ii) = dt * (dot(A3,SC)-dot(A1,SC)); % partial derivative of range w.r.t. along-track component of orbit adjustment velocity
                        orbvm(6,ii) = dt * (dot(V3,SC)-dot(V1,SC)); % partial derivative of range w.r.t. vertical    component of orbit adjustment velocity
                    else
                        orbvm(1,ii) = dot(H1,SC); % partial derivative of range w.r.t. horizontal  component of orbit adjustment vector
                        orbvm(2,ii) = dot(A1,SC); % partial derivative of range w.r.t. along-track component of orbit adjustment vector
                        orbvm(3,ii) = dot(V1,SC); % partial derivative of range w.r.t. vertical    component of orbit adjustment vector
                        orbvm(4,ii) = dt * dot(H1,SC); % partial derivative of range w.r.t. horizontal  component of orbit adjustment velocity
                        orbvm(5,ii) = dt * dot(A1,SC); % partial derivative of range w.r.t. along-track component of orbit adjustment velocity
                        orbvm(6,ii) = dt * dot(V1,SC); % partial derivative of range w.r.t. vertical    component of orbit adjustment velocity
                    end
                end
            end
            fprintf(1,'\nFinished calculating in %#10.4f seconds\n',toc(tstart));
            
            if orbopt == 1
                % Calculate orbital quantities at center for slave
                [xs, ys, zs, xdot, ydot, zdot, mjd, secs, orbnum] = readorb(ofnames2{i});
                if min(orbnum) ~= abs(islav(i))
                    warning(sprintf('Orbit numbers differ %d %d %d\n',min(orbnum),abs(islav(i)),min(orbnum)-abs(islav(i))));
                end
                [UC, DC, NC, RC, HC, AC, VC, near_mjdC, near_secC] = orbvectors2(loncenter,latcenter,0.0 ...
                    ,xs, ys, zs, xdot, ydot, zdot, mjd, secs, orbnum);
                SC = DC/norm(DC); % unit look vector target w.r.t. satellite
                for ii=i1:i2
                    if ii == i1
                        fprintf(1,'\nCalculating look vectors for SLAVE of pair %d for pixels %d to %d from orbit file %s\nPixels numbered:\n',i,i1,i2,ofnames2{i});
                        tstart=tic;
                    end
                    % calculate unit vector that varies across scene
                    if mod(ii,10) == 0; fprintf(1,'%4d ',ii); end
                    if mod(ii,100) == 0; fprintf(1,'\n'); end
                    %  calculate components of baseline
                    [U2, D2, N2, R2, H2, A2, V2, near_mjd2, near_sec2] = orbvectors2(alon(ii),alat(ii),xyzm(3,ii)...
                        ,xs, ys, zs, xdot, ydot, zdot, mjd, secs, orbnum);
                    if abs(near_mjd2 - near_mjdC) > 0
                        fprintf(1,'Longitudes %12.6f %12.6f %12.6f\n',loncenter,alon(ii),loncenter-alon(ii));
                        fprintf(1,'Latitudes  %12.6f %12.6f %12.6f\n',latcenter,alat(ii),latcenter-alat(ii));
                        dt = 86400.0*(near_mjd2 - near_mjdC) + (near_sec2 - near_secC); % time difference in seconds
                        warning(sprintf('different MJD %d %d %d\n',near_mjd2,near_mjdC,near_mjd2 - near_mjdC));
                    else
                        dt = near_sec2 - near_secC; % time difference in seconds
                    end
                    SC = D2/norm(D2); % unit look vector target w.r.t. satellite
                    
                    if ismember(pselect,[7,9]) == 1
                        %  calculate components of baseline
                        [U4, D4, N4, R4, H4, A4, V4, near_mjd4, near_sec4] = orbvectors2(alon(ii)+dlon,alat(ii),xyzm(3,ii)...
                            ,xs, ys, zs, xdot, ydot, zdot, mjd, secs, orbnum);
                        S4 = D4/norm(D4); % unit look vector target w.r.t. satellite
                        orbvs(1,ii) = dot(H4,SC) - dot(H2,SC); % partial derivative of range w.r.t. horizontal  component of orbit adjustment vector
                        orbvs(2,ii) = dot(A4,SC) - dot(A2,SC); % partial derivative of range w.r.t. along-track component of orbit adjustment vector
                        orbvs(3,ii) = dot(V4,SC)-  dot(V2,SC); % partial derivative of range w.r.t. vertical    component of orbit adjustment vector
                        orbvs(4,ii) = dt * (dot(H4,SC)-dot(H2,SC)); % partial derivative of range w.r.t. horizontal  component of orbit adjustment velocity
                        orbvs(5,ii) = dt * (dot(A4,SC)-dot(A2,SC)); % partial derivative of range w.r.t. along-track component of orbit adjustment velocity
                        orbvs(6,ii) = dt * (dot(V4,SC)-dot(V2,SC)); % partial derivative of range w.r.t. vertical    component of orbit adjustment velocity
                    else
                        orbvs(1,ii) = dot(H2,SC); % partial derivative of range w.r.t. horizontal  component of orbit adjustment vector
                        orbvs(2,ii) = dot(A2,SC); % partial derivative of range w.r.t. along-track component of orbit adjustment vector
                        orbvs(3,ii) = dot(V2,SC); % partial derivative of range w.r.t. vertical    component of orbit adjustment vector
                        orbvs(4,ii) = dt * dot(H2,SC); % partial derivative of range w.r.t. horizontal  component of orbit adjustment velocity
                        orbvs(5,ii) = dt * dot(A2,SC); % partial derivative of range w.r.t. along-track component of orbit adjustment velocity
                        orbvs(6,ii) = dt * dot(V2,SC); % partial derivative of range w.r.t. vertical    component of orbit adjustment velocity
                    end
                end
            end
            fprintf(1,'\nFinished calculating in %#10.4f seconds\n',toc(tstart));
        else
            fprintf(1,'Assuming constant unit vector [E,N,U] %10.4f %10.4f %10.4f\n'...
                ,unitv0(1),unitv0(2),unitv0(3));
            % constant across scene
            for ii=i1:i2
                unitv(1,ii)=unitv0(1);
                unitv(2,ii)=unitv0(2);
                unitv(3,ii)=unitv0(3);
            end
        end
    end
return



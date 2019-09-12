function Mw = moment_magnitude(M0inNm)
% given seismic moment M0 in Newton-meters, calculate moment magnitude Mw
% 20190813 Kurt Feigl
% 
% function Mw = moment_magnitude(M0inNm)
%
% Examples:
% moment_magnitude(1.2e9)
% ans =
%           0.02
% Chile 1960 (low estimate)
% moment_magnitude(1.4e23)
% ans =
%           9.40
% Chile 1960 (high estimate)
% moment_magnitude(2.8e23)
% ans =
%           9.60

%https://en.wikipedia.org/wiki/Moment_magnitude_scale
%The symbol for the moment magnitude scale is Mw?, with the subscript "w"
%meaning mechanical work accomplished. The moment magnitude Mw? is a
%dimensionless value defined by Hiroo Kanamori[53] as
%
%{\displaystyle M_{\mathrm {w} }={\frac {2}{3}}\log _{10}(M_{0})-10.7,} M_{\mathrm {w} }={\frac {2}{3}}\log _{10}(M_{0})-10.7,
%
% where M0? is the seismic moment in dyne?cm (10?7 N?m).[54] The constant
% values in the equation are chosen to achieve consistency with the
% magnitude values produced by earlier scales, such as the Local Magnitude
% and the Surface Wave magnitude. Thus, a magnitude zero microearthquake
% has a seismic moment of approximately 1.2×109 N?m, while the Great
% Chilean earthquake of 1960, with an estimated moment magnitude of
% 9.4?9.6, had a seismic moment between 1.4×1023 N?m and 2.8×1023 N?m.
% 

Mw = (2/3)*log10(M0inNm) - 6.03;

end


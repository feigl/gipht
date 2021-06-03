%% Prefixes
%% Using Prefixes with Unit Names
% Prefixes can be used with defined unit names to expand the scope of unit
% names. You can, for example, create kilometers from meters by adding the
% prefix "kilo".
unit(1,'kilometer')
%%
% The answer is returned as 1000 meters. The term "kilometer" is a standard
% and well recognized term. Prefixes can also be used used in conjunction
% with terms that are not as well recognized. For example, the diameter of
% the earth at the equator is
Diameter=unit(7626.28,'miles');
%%
% If the earth were perfectly spherical, the volume of the earth would be
Volume=4/3*pi*(Diameter/2)^3;
%%
% Just for fun, the volume of the earth can be expressed in
% "megateaspoons".
convert(Volume,'megateaspoons')
%% Available Prefixes and Their Values
%
%    yatto          1e24
%    zetta          1e21
%    exa            1e18
%    peta           1e15
%    tera           1e12
%    giga           1e9
%    mega           1e6
%    kilo           1e3
%    hecto          1e2
%    deka           10
%    deci           0.1
%    centi          1e-2
%    milli          1e-3
%    micro          1e-6
%    nano           1e-9
%    pico           1e-12
%    femto          1e-15
%    atto           1e-18
%    zepto          1e-21
%    yacto          1e-24
%
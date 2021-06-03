%% A script to test all functions of "unit" and "convConstant"

%%
clear all
%%
% Test that the conversion to fundamental units and back achieves the
% original result
unitType=unit.available_units;
N=numel(unitType);
pass(N)=false;
fail=false;
for I=1:N
    eval(['Unit=' unitType{I}{1},';']);
    x=convert(unit(Unit{1}),Unit{1});
    if x~=Unit{1};
        fail=true;
        disp(['The test of ' Unit{1} ' has FAILED! to convert correctly']) 
    end % if
end % for
if ~fail
    disp('     All circular conversion tests PASS!')
end % if
%%
% Check addition and subtraction of unit
sum=-unit('yd')+unit(2,'yds');
diff=sum-unit('yd');
if diff==unit(0,'yds')
    disp('     Addition and subtraction PASS!')
else
    disp('     Addition and subtraction FAIL!')
end % if
%%
% Check addition and subtraction of dimensionless units
x=unit(3,'');
if 2*x-3+x==6
    disp('     Dimensionless addition and subtraction PASS!')
else
    disp('     Dimensionless addition and subtraction FAIL!')
end % if  

%%
% Check division and equality
if (unit('yd')/unit('yds'))==unit('ft')/unit('ft')
    disp('     Division and equality PASS!')
else
    disp('     Division and equality FAIL!')
end
%%
% Check multiplication and exponentation
if unit('yd')*unit('yd*yd')==unit('yd*yd^3*yd')/unit('yd^2')
    disp('     Multiplication and exponentiation PASS!')
else
    disp('     Multiplication and exponentiation FAIL!')
end % if
%%
% Check reciprocals
if 1/unit('yd')==unit('1/yd')
    disp('     Reciprocals PASS!')
else
    disp('     Reciprocals FAIL!')
end
%%
% Test for equality of empty units
if unit('yd')/unit('yd')==unit('')
    disp('     Empty unit test PASS!')
else
    disp('     Empty unit test FAIL!')
end % if
%% 
% Check sin and cos
if sin(unit('yds')/unit(2,'yds'))^2+cos(unit('yds')/unit(2,'yds'))^2==1
    disp('     Sin and cos functions PASS!')
else
    disp('     Sin and cos functions PASS!')
end % if
%%
% Check square root
if sqrt(unit('yds^4'))==unit('yds^2')
    disp('     Square root PASSES!')
else
    disp('     Square root FAILS!')
end % if
%%
% Check raising a unit to a array power
if unit('yd').^2==unit('yd^2') %#ok<BDSCA>
    disp('     Raising a Unit to a power PASSES')
else
    disp('     Raising a Unit to a power FAILS')
end
%% 
% Check sin and csc
if csc(unit('yds')/unit(2,'yds'))==sin(unit('yds')/unit(2,'yds'))^-1
    disp('     Sin and csc functions PASS!')
else
    disp('     Sin and csc functions PASS!')
end % if
%%
% Check cos and sec
if sec(unit('yds')/unit(2,'yds'))==cos(unit('yds')/unit(2,'yds'))^-1
    disp('     Cos and sec functions PASS!')
else
    disp('     Cos and sec functions PASS!')
end % if
%%
% Check tangent and cotangent
if tan(unit('yds')/unit(2,'yds'))*cot(unit('yds')/unit(2,'yds'))==1
    disp('     Tan and cot functions PASS!')
else
    disp('     Tan and cot functions FAIL!')
end % if
%%
% Test radians with sin and cos
x=unit(pi/4,'radians');
y=convert(x,'radians');
if sin(y)^2+cos(y)^2==1
    disp('     Radians with sin and cos PASS!')
else
    disp('     Radians with sin and cos FAIL!')
end
%%
% Check the inverse trig functions
x=sin(asin(unit))+cos(acos(unit))+tan(atan(unit))+cot(acot(unit))+sec(asec(unit))+csc(acsc(unit));
if x==6
    disp('     Inverse trig functions PASS!')
else
    disp('     Inverse trig functions FAIL!')
end % if
%%
% Check hyperbolic sin and cos
if cosh(unit(pi/2,'rad'))==sqrt(1+sinh(unit(pi/2,'rad'))^2)
    disp('     Hyperbolic sin and cos PASS!')
else
    disp('     Hyperbolic sin and cos FAIL!')
end
%%
% Check hyperbolic tan and cot
if coth(unit)==1/tanh(unit)
    disp('     hyperbolic tan and cot PASS!')
else
    disp('     hyperbolic tan and cot FAIL!')
end
%%
% Check hyperbolic se and cosec
if abs(csch(unit)-sech(unit)/sqrt(1-sech(unit)^2))<1e-15
    disp('     Hyperbolic sec and csc PASS!')
else
    disp('     Hyperbolic sec and csc FAIL!')
end % if
%%
% Check the inverse hyperbolic functions
x=sinh(asinh(unit))+cosh(acosh(unit))+tanh(atanh(unit))+coth(acoth(unit))+sech(asech(unit))+csch(acsch(unit));
if x==6
    disp('     Inverse hyperbolic functions PASS!')
else
    disp('     Inverse hyperbolic functions FAIL!')
end % if
%%
% Check conversion from type "double"
if unit(3)==3*unit('')
    disp('     Conversion from type "double" PASS!')
else
    disp('     Conversion from type "double" FAIL!')
end
%%
% Check natural log and exponential
if log(exp(unit(3,'yds')))-exp(log(unit(3,'yds')))<1e-15
    disp('     Natural log and exponential PASS!')
else
    disp('     Natural log and exponential FAIL!')
end
%%
% Check common log and power of ten
x=unit(3,'yds');
if 10^log10(x)-x.value<1e-15
    disp('     Common log and power of 10 PASS!')
else
    disp('     Common log and power of 10 FAIL!')
end % if
%%
% Check that the fundamental units all execute correctly
unit('m*s*kelvin*candela*kg^3*mole*amp^2');
unit('ft^2');
%%
% Check that all matrix functions execute without returning an error
% (validity of the result is not checked)
x=magic(5)*unit('newtons');
x';x.';-x;+x;x+x;x-x;x*x;x.*x;x./x;x.\x;x/x;x\x; %#ok<MNEFF,VUNUS>
x==x;x~=x;x>x;x<x;x<=x;x>=x; %#ok<VUNUS,EQEFF>
disp('     Matrix operations execute without returning an error!')
%%
% Test the display routine (Check visually)
unit(10,'yds')
unit(1:10,'yds')
%% 
% Test the search function (Check visually)
unit.search('farad')
%% Stream flow example
clear all
%% Load the date time and flow data from the three files
fid=fopen('glenwood.txt');
GlenFlow=textscan(fid,'%*s %*s %s %s %*s %*s %f','CollectOutput',true);
fclose(fid);
fid=fopen('kremmling.txt');
Kremflow=textscan(fid,'%*s %*s %s %s %*s %*s %f','CollectOutput',true);
fclose(fid);
fid=fopen('dotsero.txt');
Dotsflow=textscan(fid,'%*s %*s %s %s %*s %*s %f','CollectOutput',true);
fclose(fid);
%% Create a date number from the date and time. Use the Dotsero data for the set as their is one less sample in this set
[samples,~]=size(Dotsflow{1});
dateTime{samples}='';
% Concatinate the the two cells holding date and time
for N=1:samples
    dateTime{N}=[Dotsflow{1,1}{N,1} ' ' Dotsflow{1,1}{N,2}];
end
D=datenum(dateTime);
%% Convert the date numbers to units of days
Date=convert(unit(D,'day'),'days');
%% Create the flow information with units for plotting
Flow.Dotsero=unit(Dotsflow{1,2}(1:N),'ft^3/sec');
Flow.Kremmling=unit(Kremflow{1,2}(1:N),'ft^3/sec');
Flow.Glenwood=unit(GlenFlow{1,2}(1:N),'ft^3/sec');
%% Plot the flow in metric units of cubic meters per second
% Return the legend handle so the legend can be positioned to not obscure
% data
[~,legHandle]=plot(Date,Flow);
%% Make the x axis display dates as month/day
datetick('x','mm/dd')
%% Title the graph
title('2010 Flow Rate for the Colorado River at Three Gage Stations')
%% Position the legend so it does not obscure data
% Data obtained by positioning the legend and then getting the position
% information manually.
set(legHandle,'Position',[0.6568 0.5885 0.2143 0.1444]);
%% Convert the data to acre feet per day and replot 
Flow.Dotsero=convert(Flow.Dotsero,'acre*feet/day');
Flow.Kremmling=convert(Flow.Kremmling,'acre*feet/day');
Flow.Glenwood=convert(Flow.Glenwood,'acre*feet/day');
%% Plot the flow in acre-feet per day
[~,legHandle]=plot(Date,Flow);
%% Make the x axis display dates as month/day
datetick('x','mm/dd')
%% Title the graph
title('2010 Flow Rate for the Colorado River at Three Gage Stations')
%% Position the legend so it does not obscure data
set(legHandle,'Position',[0.6568 0.5885 0.2143 0.1444]);
%% Plot a single vector of x and y units
Dotsero=Flow.Dotsero;
plot(Date,Dotsero)
datetick('x','mm/dd')
title('2010 Flow Rate of the Colorado River at Gage Station Dotsero')
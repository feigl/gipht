classdef convConstant
    % Maintains the conversion constants used by the class "unit"
    %
    %   An excellent reference, by Prof. Russ Rowlett, University of
    %   North Carolina at Chapel Hill, is available on the Web at:
    %   http://www.unc.edu/~rowlett/units
    %
    % AUTHOR:
    %           John McDermid
    % CREATED:
    %           December, 2010
    % REVISED:
    %           January 20, 2012   Correct the definition of weber
    
    properties
        unitName       % The name of the unit
        convValue      % When the unit is converted to the conversion unit, this is the conversion value
        convUnit       % The name of the unit being converted to
        fundamental    % A boolean representing whether this a fundamental unit
    end % properties
    
    methods
        function obj=convConstant
            % The convConstant constructer method
            % Syntax:
            %       obj=convConstant
            %           where obj.unitName='none'
            %                 obj.convValue=[1 0]
            %                 obj.convUnit=none
            %                 obj.fundamental is empty
            
            obj.unitName='none';obj.convValue=[1 0];obj.convUnit='none';obj.fundamental='';  
        end % function
        
        function y=define(obj,x)
            % This method defines a unit in terms of its conversion value
            % and conversion unit.
            % Syntax:
            %       y=define(obj,x)
            %           where x is the string to be defined
            
            y=obj;
            switch x
                % Fundamental Units (International Units or SI)
                case {'seconds','second','sec','s'}          % time in seconds
                    y.unitName='s';y.convValue=[1 0];y.convUnit='s';y.fundamental=true;
                case {'meter','meters','m'}                  % length in meters
                    y.unitName='m';y.convValue=[1 0];y.convUnit='m';y.fundamental=true; 
                case {'kg'}                                  % mass in kilograms
                    y.unitName='kg';y.convValue=[1 0];y.convUnit='kg';y.fundamental=true;
                case {'ampere','amperes','amp','amps','A'}   % electrical current in amperes
                    y.unitName='A';y.convValue=[1 0];y.convUnit='A';y.fundamental=true;
                case {'kelvin','K'}                          % temperature in degrees kelvin
                    y.unitName='K';y.convValue=[1 0];y.convUnit='K';y.fundamental=true;
                case {'mole','moles','mol'}                  % amount of a substance in moles
                    y.unitName='mol';y.convValue=[1 0];y.convUnit='mol';y.fundamental=true;
                case {'candela','candelas','cd'}             % light intensity in candelas
                    y.unitName='cd';y.convValue=[1 0];y.convUnit='cd';y.fundamental=true;
                
                % Derived Units
                case {'newton','newtons','N'}                % force in newtons
                    y=units2convFac(unit('kg*m/s^2'));
                case {'pascal','pascals','Pa'}               % pressure in pascals
                    y=units2convFac(unit('N/m^2'));
                case {'joule','joules','J'}                  % energy in joules
                    y=units2convFac(unit('N*m'));
                case {'watt','watts','W'}                    % power in watts
                    y=units2convFac(unit('J/s'));
                case {'coulomb','coulombs','C'}              % electrical charge in coulombs
                    y=units2convFac(unit('A*s'));
                case {'volt','volts','V'}                    % voltage in volts
                    y=units2convFac(unit('joule/C'));
                case {'weber','webers'}                      % magnetic flux in webers
                    y=units2convFac(unit('volt*second'));
                case {'tesla','teslas'}                      % magnetic flux density in teslas
                    y=units2convFac(unit('weber/m^2'));
                case {'gauss'}                               % magnetic flux density in gauss
                    y=units2convFac(1e-4*unit('tesla'));
                case {'farad','farads','F'}                  % capacitance in farads
                    y=units2convFac(unit('coulombs/volt'));
                case {'ohm','ohms'}                          % resistance in ohms
                    y=units2convFac(unit('volt/amp'));
                case {'henry','henrys','H'}                  % inductance in henrys
                    y=units2convFac(unit('volts*seconds/amps'));
                case {'diopter','diopters'}                  % focal length in diopters
                    y=units2convFac(1/unit('m'));
                    
                % CGS
                case {'cm'}                                  % length in centimeters
                    y=units2convFac(unit('centimeter'));
                case {'gram','grams','gr'}                   % weight in grams
                    y=units2convFac(unit(1/1000,'kg'));
                case {'calorie','calories','cal'}            % energy in calories
                    y=units2convFac(unit(4.1868,'joules'));
                case {'Calorie','Calories','kcal'}           % energy in kilocalories
                    y.unitName='kcal';y.convValue=[1000 0];y.convUnit='cal';y.fundamental=false;
                case {'dyne','dynes'}                        % energy in dynes
                    y=units2convFac(unit(1e-5,'newtons'));
                    
                % Energy
                case {'Btu','BTU'}                           % energy in british thermal units
                    y=units2convFac(1055.056*unit('J'));
                case {'eV'}                                  % energy in electron volts
                    y=units2convFac(160.217646e-12*unit('J'));
                    
                % Power
                case {'hp'}                                  % power in horse power
                    y=units2convFac(745.699*unit('watt'));
                    
                % Length
                case {'km'}                                  % length in kilometers
                    y.unitName='km';y.convValue=[1000 0];y.convUnit='m';y.fundamental=false;
                case {'yard','yards','yd','yds'}             % length in yards
                    y.unitName='yd';y.convValue=[3 0];y.convUnit='feet';y.fundamental=false;
                case {'foot','feet','ft'}                    % length in feet
                    y.unitName='ft';y.convValue=[12 0];y.convUnit='inch';y.fundamental=false;
                case {'inch','inches','in'}                  % length in inches
                    y.unitName='in';y.convValue=[0.0254 0];y.convUnit='meter';y.fundamental=false;
                case {'mile','miles'}                        % length in miles
                    y.unitName='mile';y.convValue=[5280 0];y.convUnit='feet';y.fundamental=false;
                case {'fathom','fathoms'}                    % length in fathoms
                    y.unitName='fathom';y.convValue=[6 0];y.convUnit='feet';y.fundamental=false;
                case {'rod','rods'}                          % length in rods
                    y.unitName='rod';y.convValue=[5.5 0];y.convUnit='yards';y.fundamental=false;
                case {'furlong','furlongs'}                  % length in furlongs
                    y.unitName='rod';y.convValue=[40 0];y.convUnit='rod';y.fundamental=false;
                case {'nmi','naut_mi'}                       % length in nautical miles
                    y.unitName='nmi';y.convValue=[6076.11549 0];y.convUnit='feet';y.fundamental=false;
                case {'league','leagues'}                    % length in leagues
                    y.unitName='league';y.convValue=[3 0];y.convUnit='nmi';y.fundamental=false;
                case {'angstrom'}                            % length in angstroms
                    y.unitName='angstrom';y.convValue=[1e-10 0];y.convUnit='meter';y.fundamental=false;
                case {'lightyear'}                           % length in light years
                    y=units2convFac(unit('c*yr'));
                case {'au'}                                  % length in astronomical units
                    y=units2convFac(149597870*unit('kilometer'));
                case {'parsec','parsecs'}                    % distance that a star appears to shift when the earth moves one "au"
                    y=units2convFac(206264.8*unit('au'));
                    
                % Time
                case {'minute','minutes','min'}              % time in minutes
                    y.unitName='min';y.convValue=[60 0];y.convUnit='s';y.fundamental=false;
                case {'hour','hours','hr'}                   % time in hours
                    y.unitName='hr';y.convValue=[60 0];y.convUnit='min';y.fundamental=false;
                case {'day','days'}                          % time in days
                    y.unitName='day';y.convValue=[24 0];y.convUnit='hr';y.fundamental=false;
                case {'week','weeks'}                        % time in weeks
                    y.unitName='week';y.convValue=[7 0];y.convUnit='day';y.fundamental=false;
                case {'fortnight','fortnights'}              % time in fortnights
                    y=units2convFac(unit(2,'weeks'));
                case {'year','years','yr'}                   % time in years
                    y=units2convFac((365+97/400)*unit('day')-26.0237*unit('sec'));
                case {'century','centuries'}                 % time in centuries
                    y=units2convFac(unit(100,'yr'));
                    
                % Mass
                case {'pound','pounds','lb','lbs','lbm'}     % weight in pounds (avoirdupois)
                    y.unitName='lb';y.convValue=[0.45359237 0];y.convUnit='kg';y.fundamental=false;
                case {'ounce','ounces','oz'}                 % weight in ounces (avoirdupois)
                    y.unitName='oz';y.convValue=[1/16 0];y.convUnit='lb';y.fundamental=false;
                case {'grain'}                               % weight in grains
                    y.unitName='grain';y.convValue=[1/7000 0];y.convUnit='lb';y.fundamental=false;
                case {'ton'}                                 % weight in tons
                    y.unitName='ton';y.convValue=[2000 0];y.convUnit='lb';y.fundamental=false;
                case {'amu','u'}                             % atomic mass units
                    y=units2convFac(unit('gram')/unit('Avogadro*mole'));
                    
                % Force
                case {'lbf'}                                 % pounds force
                    y=units2convFac(unit('lb*g'));
                
                % Temperature
                case {'celsius','degC'}                      % temperature in degrees Celsius
                    y.unitName='degC';y.convValue=[1 273.15];y.convUnit='kelvin';y.fundamental=false;
                case {'degF'}                                % temperature in degrees Fahrenheit
                    y.unitName='degF';y.convValue=[5/9 -32*5/9];y.convUnit='degC';y.fundamental=false;
                case {'degR'}                                % temperature in degrees Rankine
                    y.unitName='degR';y.convValue=[5/9 0];y.convUnit='kelvin';y.fundamental=false;
                    
                % Area
                case {'acre','acres'}                         % area in acres
                    y=units2convFac(160*unit('rod^2'));
                case {'are'}                                  % area in "are" (pronounced "air")
                    y=units2convFac(100*unit('m^2'));
                case {'hectare','hectares'}                   % area in hectares
                    y=units2convFac(100*unit('are'));
                case {'barn','barns'}                         % area in barns
                    y=units2convFac(1.e-28*unit('m^2'));
                    
                % Volume
                case {'liter','liters','L'}                   % volume in liters
                    y=units2convFac(unit('m^3')/1000);
                case {'cc'}                                   % volume in cubic centimeters
                    y=units2convFac(unit('L')/1000);
                case {'fldoz'}                                % volume in fluid ounces
                    y=units2convFac(0.029573531*unit('L'));
                case {'tablespoon','tablespoons','tblsp'}     % volume in table spoons
                    y=units2convFac(unit('fldoz')/2);
                case {'teaspoon','teaspoons','tsp'}            % volume in tea spoons
                    y=units2convFac(unit('tblsp')/3);
                case {'cup','cups'}                           % volume in cups
                    y=units2convFac(8*unit('fldoz'));
                case {'pint','pints'}                         % volume in pints
                    y=units2convFac(2*unit('cup'));
                case {'quart','quarts'}                       % volume in quarts
                    y=units2convFac(2*unit('pint'));
                case {'gallon','gallons'}                     % volume in gallons (U.S.)
                    y=units2convFac(4*unit('quart'));
                case {'acre_foot','acre_feet'}                % volume in acre feet
                    y=units2convFac(unit('acre*feet'));
                    
                % Velocity
                case {'mph'}                                  % velocity in miles per hour
                    y=units2convFac(unit('miles/hour'));
                case {'mach'}                                 % velocity in mach numbers (speed of sound)
                    y=units2convFac(741.8*unit('miles/hour'));
                case {'knot','knots','kn'}                    % velocity in nautical miles per hour
                    y=units2convFac(unit('nmi/hour'));
                    
                    
                % Physical Constants
                case {'c'}                                    % the speed of light in a vacuum
                    y=units2convFac(unit(299792458,'m/s'));
                case {'g'}                                    % acceleration due to gravity
                    y=units2convFac(unit(9.80665,'meter/second^2'));
                case {'e'}                                    % charge on an electron
                    y=units2convFac(unit(1.6021773349e-19,'C'));
                case {'me'}                                   % mass of an electron
                    y=units2convFac(unit(9.109389754e-31,'kg'));
                case {'mp'}                                   % mass of a proton
                    y=units2convFac(unit(1.672623110e-27,'kg'));
                case {'Planck','h'}                           % energy expended over time
                    y=units2convFac(unit(6.62607540e-34,'joule*sec'));
                case {'alpha'}                                % fine structure constant
                    y=units2convFac(2*pi*1.e-7*unit('c*e^2/h'));
                case {'R'}                                    % gas constant
                    y=units2convFac(8.31451070*unit('J/mol/K'));
                case {'Avogadro','NA'}                        % Avogadro's number
                    y=units2convFac(6.022136736e23/unit('mol'));
                case {'muZero'}                               % magnetic permeability of free space
                    y=units2convFac(unit(4*pi*1e-7,'N*A^-2'));
                case {'epsilonZero'}                          % dielectric constant of free space
                    y=units2convFac(unit(8.854187817e-12,'F/meter'));
                    
                % Pressure
                case {'atm'}                                  % pressure in atmospheres
                    y=units2convFac(101325*unit('Pa'));
                case {'torr'}                                 % pressure in torr
                    y=units2convFac(133.322*unit('pascal'));
                case {'bar','bars'}                           % pressure in bars
                    y=units2convFac(100*unit('kilopascals'));
                    
          % Other
                case {'radian','radians','rad'}              % angle in radians
                    y=units2convFac(unit);  
                case {'degree','degrees','deg'}              % angle in degrees
                    y=units2convFac(pi/180*unit('rad'));
                case {'cycle','cycles','rev'}                % a periodic unit
                    y=units2convFac(2*pi*unit('rad')); 
                case {'hertz','Hz'}                          % frequency in hertz
                    y=units2convFac(unit('cycles/second'));
                    
                otherwise
                    error(['Unit known as "' x '" is not in the conversion table!'])
            end % switch
        end % function
        
        function newObj=conv2fund(obj)
            % Converts an object to an object with fundamental units
            % Syntax:
            %         newObj=conv2fund(obj)
            
            if obj.fundamental==false
                X=define(convConstant,obj.convUnit);
                Y=conv2fund(X);
                newObj=Y;
                newObj.convValue=combineValues(Y.convValue,obj.convValue);
            else
                newObj=obj;
            end % if
        end % function
        
    end % methods
    
end % classdef

function newVal=combineValues(value,convValue)
% If Y1=A*X1+B and Y2=C*X2+D and Y1 is substituted for the value of X2,
% then Y2=C*(A*X1+B)+D=(C*A*X1+C*B)+D. This operation is needed where a
% conversion of units also has an offset value, such as Centigrade to
% Fahrenheit. 
C=value(1);D=value(2);
A=convValue(1);B=convValue(2);
% The value of C*A is determined.
a=C*A;
% The value of (C*B+D) is determined.
b=C*B+D;
% The new values are assigned to the output array.
newVal=[a b];
end

function y=units2convFac(x)
% Converts a "unit" object to a conversion constant
if isa(x,'unit')
    y=convConstant;
    y.unitName=x.name;
    y.convValue=[x.value x.offset];
    y.convUnit=x.name;
    y.fundamental=true;
else
    error('Argument "x" of "units2convFac(x)" must be of type "unit"')
end % if
end % function
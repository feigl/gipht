%% Defining Variables and their Units
%% Fundamental Units
% A "unit" statement accepts a number and a unit name and converts these
% to international units (SI). There are just seven different units in the
% international system of units. These are:
%
%           seconds (s)    -- units of time
%           meters (m)     -- units of length
%           kilograms (kg) -- units of mass
%           ampere (A)     -- units of electrical current
%           kelvin (K)     -- units of temperature
%           mole (mol)     -- units of substance quantity
%           candela (cd)   -- units of light intensity
%
% The parenthesis contain the abreviation for the fundamental unit.
% Abbreviations are used as the unit names. These names may appear as a
% simple abbreviation such as when specifying a yard.
unit(1,'yard')
%%
% or as a mathematical combination such as when specifying an ohm.
unit(1,'ohm')
%%
%%
%% Finding Defined Unit Names
% The available units and their names can be found with the "unit.search"
% command. To see if there is a unit for pounds, search for "pound" with
% the search command.
unit.search('pound')
%%
% Two choices are presented, one for pounds avoirdupois and one for pounds
% force. The search command returns all descriptions that contain "pound".
% If you would like to see the entire table of constants, you can use the
% word "all".
unit.search('all')
%%
%%
%% Creating New Units from Existing Unit Names
% The existing unit names will never cover all of the instances of unit
% names you may require. For example, the conversion table contains no name
% for the velocity in feet per second. However, it is possible to construct
% a unit name through mathematical operators. The ballistic table for the
% 0.17 Remington Fireball cartridge shows the velocity to be
Velocity=unit([4250 3594 3028 2529 2081],'feet/second');
%%
% where the velocity is expressed in 100 yard increments from the muzzle
% out to 400 yards. The velocity unit name is defined as the ratio of units
% "feet" and "second".
%%
% The operands available to combine defined units are *, /, and ^. Addition
% and subtraction are not permitted with unit names. The number 1 is
% permitted so that reciprocals can be formed. For example
unit('1/yd')
%%
%%
%% Converting to an Equivalent Unit of Measure
% The "unit" statement always expresses the result in the SI system of
% measure. Units can be converted to any equivalent unit of measure with
% the "convert" statement. The velocity, defined above, can be expressed in
% mach numbers.
convert(Velocity,'mach')
%%
% The "unit" statement reports unit names with the shortest abreviation.
% Unit names expressed in a convert statement are reported exactly the way
% they are expressed. For example
x=unit(1,'yds')
x=convert(x,'meters')
%%
% The unit names in a "convert" statement can also be created with
% operators. The velocity in the example above can be expressed in
% kilometers per hour by
convert(Velocity,'kilometers/hour')
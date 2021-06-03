%% Adding Your Own Conversion Constants
%% Where to Add
% The conversion constants are contained within the "define" method in the
% "convConstant" object. Each constant is part of a series of case
% statements from one large switch statement. Here is an example case
% statement.
%
%  case {'yard','yards','yd','yds'}             % length in yards
%    y.unitName='yd';y.convValue=[3 0];y.convUnit='feet';y.fundamental=false;
%
% New conversion constants are added by duplicating the structure above.
%% Syntax
% The strings within the curly brackets are the permissible names of the
% units. The convention is to add the full name in both singular and plural
% forms and abbreviations for the units if they exist.
%
% The comment is used as the search material for the statement
% "unit.search". It should start with a lower case letter and complete the
% statement "Units of ...".
% 
% The next line contains four statements. The property "y.unit" contains
% the unit name. The property "y.convValue" is the pair of constants that
% are applied to the conversion unit to get the unit name. For example,
%
%  yard=[3 0]*feet
% 
% The first number multiplies the conversion unit and the second adds an
% offset if required. Most of the conversion constants have an offset of
% zero. An exception is Fahrenheit to Celsius conversion.
%
% The conversion unit "y.convUnit" is the name of the conversion unit. In
% the example above it is "feet". This can only be a simple name (no
% operaters allowed). It must be a unit name that already exists in the
% define method.
%
% The statement "y.fundamental" is a boolean variable that describes
% whether this is a member of the SI system of units. Since all members of
% the fundamental SI units have already been described, this value is
% always false.
%
% Sometimes it is necessary to use mathematical operators to construct a
% new unit. For example, to define a light year.
%
%  case {'lightyear'}                           % length in light years
%    y=units2convFac(unit('c*yr'));
%
% The "unit" constructor is used to multiply the speed of light "c" by the
% time in years. The function "units2convFac" is then used to directly
% create the structure for y.
%% Procedure
% First, demonstrate that conversion unit already exists in the define
% method. As an example, prove that "feet" exists by typing
unit('feet')
%%
% This should return an answer in meters. Next, demonstrate that the unit
% you wish to define does not already exist as the name of some other unit.
% Suppose that you wish to define a unit called "shortmile". Type
unit('shortmile')
%%
% An error should be reported stating that the unit is not in the
% conversion table. Check each equivalent unit name to be sure that no
% conflict exists.
%
% Finally, add the case and definition statements to the "define" method
% and save the file. Then type "clear classes" and type "unit('<new
% unit>'). The value should be returned in fundamental units. Use the
% "unit.search" statement and search for a word that is in your comment.
% One of the pieces of information returned should be your new unit.

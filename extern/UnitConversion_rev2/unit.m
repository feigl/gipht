classdef unit
    % "unit" attaches physical units to a variable.
    % Examples:
    %       unit(1,'mile')
    %
    %       ans =
    %             1609.34 m
    %
    %       y=unit('furlongs/fortnight')
    %
    %       ans =
    %             0.00016631 m/s
    %
    % Any variable defined by "unit" can be converted to to a compatible
    % unit (you cannot convert length to area for example) by using the
    % "convert" method. The next example converts the variable y, defined
    % above, to feet per minute.
    %       convert(y,'feet/minute')
    %
    %       ans =
    %             0.0327381 feet/minute
    % Syntax:
            %   U=unit        -- called with no parameters it returns a unit
            %                     value of one and an empty name string
            %   U=unit(x)     -- Where x is of type "double" returns a unit
            %                     value of x and an empty name string
            %   U=unit('x')   -- where x is a valid unit string returns a unit
            %                     value of 1 and the unit string converted to
            %                     fundamental units
            %   U=unit(x)     -- where x is a class of "unit' simply returns x
            %                     (in fundamental units if it had been converted
            %                      to other unit values)
            %   U=unit(n,'x') -- where n is a double and x is a valid unit
            %                     string returns unit value of n and a unit
            %                     string that is x converted to fundamental
            %                     units
            %
            % Fundamental units are SI units. There are seven fundamental
            % units in the SI system of units. They are:
            %   seconds   -- units of time
            %   meters    -- units of length
            %   kilograms -- units of mass
            %   ampere    -- units of electrical current
            %   kelvin    -- units of temperature
            %   mole      -- units of substance quantity
            %   candela   -- units of light intensity
            %
            % All mathematical operations are performed in fundamental
            % units. Units can be converted to whatever units are desired
            % by using the "convert" method.
            % Syntax:
            %       U=convert(A,'B')
            %                       where A is any variable of class
            %                       "unit" and B is the desired unit. 
            %
            % The value of B may be selected by using the unit.search(str)
            % function or by combining units with operators such as the
            % string 'feet/minute'.
            %
            % Units strings can be found using unit.search('string') where
            % the string is any string that might describe the unit. Units
            % strings can also be created by using operators within the
            % unit string. There is no unit string that represents feet per
            % second like there is for mile per hour (mph). However, a unit
            % string can be constructed by using 'feet/second' as a unit
            % string. The available operators inside a unit string are *,
            % /, and ^.
            %
            % Variables created by "unit" can be make use of the following
            % mathematical operators:
            %       Unary plus, Unary minus, +,-,*,.*,/,\,./,.\,^,.^,',.'
            %       and sqrt
            % The following triginometric function will accept variables
            % with units:
            %       sin, cos, tan, cot, sec, cosec, asin, acos, atan, acot,
            %       asec, and acosec
            % The following hyperbolic functions are also available:
            %       sinh, cosh, tanh, coth, sech, cosech, asinh, acosh,
            %       atanh, acoth, acosech, and asech
            % Exponential functions are also available. They are:
            %       log10, log, and exp
            % Logical functions are also available to compare variables
            % with units. They are:
            %       ==, ~=, >, <, >=, and <=
            % There are also four plotting functions that accept variables
            % with units. They are:
            %       plot, semilogx, semilogy, and loglog
%
% AUTHOR:
%        John McDermid
% CREATED:
%        December, 2010
% REVISED:
%        January 20,2012    Correct an error in "plus" and "minus" for
%        dimensionless units
    
    properties
        value             % value represented by the unit
        offset            % offset required by the unit
        name              % string containing the unit description
        hasBeenConverted  % boolean value describing whether it was generated with the "convert" statement
    end % properties
    
    methods
        
        function U=unit(varargin)
            % The constructor function for class "unit"
            % Syntax:
            %   U=unit        -- called with no parameters it returns a unit
            %                     value of one and an empty name string
            %   U=unit(x)     -- Where x is of type "double" returns a unit
            %                     value of x and an empty name string
            %   U=unit('x')   -- where x is a valid unit string returns a unit
            %                     value of 1 and the unit string converted to
            %                     fundamental units
            %   U=unit(x)     -- where x is a class of "unit' simple returns x
            %                     (in fundamental units if it had been converted
            %                      to other unit values)
            %   U=unit(n,'x') -- where n is a double and x is a valid unit
            %                     string returns unit value of n and a unit
            %                     string that is x converted to fundamental
            %                     units
            %
            % Fundamental units are SI units. There are seven fundamental
            % units in the SI system of units. They are:
            %   seconds   -- units of time
            %   meters    -- units of length
            %   kilograms -- units of mass
            %   ampere    -- units of electrical current
            %   kelvin    -- units of temperature
            %   mole      -- units of substance quantity
            %   candela   -- units of light intensity
            %
            % All mathematical operations are performed on fundamental
            % units. Units can be converted to whatever units are desired
            % by using the "convert" method.
            %
            % Units strings can be found using unit.search('string') where
            % the string is any string that might describe the unit. Units
            % strings can also be created by using operators within the
            % unit string. There is no unit string that represents feet per
            % second like there is for mile per hour (mph). However, a unit
            % string can be constructed by 'feet/second' as a unit string.
            % The available operators inside a unit string are *, /, and ^.
            
            switch nargin
                % If no arguments are provided, return a value of unity
                % with no units attached
                case 0
                    U.value=1;
                    U.offset=0;
                    U.name='';
                    U.hasBeenConverted=false;
                case 1
                    % If there is no Unit string, create and empty unit
                    if isempty(varargin{1})
                        U.value=1;
                        U.offset=0;
                        U.name='';
                        U.hasBeenConverted=false;
                        % If there is just one argument and it is already a
                        % Unit, just return it. If its not in fundamental units
                        % then change it to fundamental
                    elseif isa(varargin{1},'unit')
                        U=varargin{1};
                        if U.hasBeenConverted
                            U=unit(U.value+U.offset,U.name);
                            U.hasBeenConverted=false;
                        else
                            U=varargin{1};
                        end % if
                        % If there is just one argument and it is a number
                        % (signed or unsigned) convert it to a number, assign
                        % the number to the value, and attach empty units.
                    elseif ~isnan(str2double(varargin{1}))
                        U.value=str2double(varargin{1});
                        U.offset=0;
                        U.name='';
                        U.hasBeenConverted=false;
                        % If varargin{1} is just a number create a Unit
                        % with an empty units string
                    elseif isa(varargin{1},'double')
                        U.value=varargin{1};
                        U.offset=0;
                        U.name='';
                        U.hasBeenConverted=false;
                    elseif isa(varargin{1},'char')
                        % If the argument is a string of characters, assign the
                        % input argument to a string and find the first token
                        % (separaters are the operators), the operator, and the
                        % remaining characters.
                        remainStr=varargin{1};
                        I=1;
                        operator='';
                        % Find all the tokens and the operators in the string (remainStr) and
                        % provide the count of operators and tokens
                        while ~isempty(remainStr)
                            [token{I},remainStr]=strtok(remainStr,'*/^');  %#ok<STTOK>
                            if ~isempty(remainStr)
                                operator{I}=remainStr(1);
                            end
                            I=I+1;
                        end
                        Ntok=numel(token);
                        Nop=numel(operator);
                        if Ntok==1
                            T1=token{1};
                            [preVal,T1]=removePrefix(T1);
                            C=define(convConstant,T1);
                            C=conv2fund(C);
                            U.value=preVal*C.convValue(1);
                            U.offset=C.convValue(2);
                            U.name=C.unitName;
                        else
                            % Initialize the variables required to build strings that can be evaluated
                            % using the "eval" statement
                            S=token{1};
                            I=1;
                            % Build and evaluate the necessary  "unit" statements
                            while I<Ntok
                                if I<=Nop-1
                                    % Look ahead to the operator after the next operator. If it is an
                                    % exponent symbol(^), perform this operation first and then apply
                                    % the next operator.
                                    if char(operator(I+1))=='^'
                                        if isa(S,'unit')
                                            S=eval(['unit(S)' operator{I} ['unit(''' token{I+1} ''')' operator{I+1} 'unit(''' token{I+2} ''')']]);
                                        else
                                            S=eval(['unit(''' S ''')' operator{I} ['unit(''' token{I+1} ''')' operator{I+1} 'unit(''' token{I+2} ''')']]);
                                        end %if
                                        I=I+1;
                                    else
                                        if isa(S,'unit')
                                            S=eval(['unit(S)' operator{I} ['unit(''' token{I+1} ''')']]);
                                        else
                                            S=eval([['unit(''' S ''')'] operator{I} ['unit(''' token{I+1} ''')']]);
                                        end %if
                                    end %if
                                elseif I==Nop
                                    % If the last operator is an exponentiation symbol(^), perform the
                                    % exponentiation operation and combine with the previous string
                                    if char(operator(Nop))=='^'
                                        if isa(S,'unit')
                                            S=eval(['unit(S)' operator{I} ['unit(' token{I+1} ')']]);
                                        else
                                            S=eval([['unit(''' S ''')'] operator{I} ['unit(''' token{I+1} ''')']]);
                                        end %if
                                    else
                                        if isa(S,'unit')
                                            S=eval(['unit(S)' operator{I} ['unit(''' token{I+1} ''')']]);
                                        else
                                            S=eval([['unit(''' S ''')'] operator{I} ['unit(''' token{I+1} ''')']]);
                                        end %if
                                    end % if
                                end %if
                                I=I+1;
                            end % while
                            U=S;
                        end % if
                        % Ensure that each unit string is expressed in the
                        % simplest form.
                        x=simplifyFund(U.name);
                        y=fundUnits2str(x);
                        U.name=y;
                        U.hasBeenConverted=false;
                    else
                        error('Parameter must be a string!');
                    end % if
                case 2
                    % If there are two variables, the first a number and
                    % the second a string, express the result as the number
                    % times the Unit.
                    if isa(varargin{1},'double') && isa(varargin{2},'char')
                        X=unit(varargin{2});
                        U=X;
                        U.value=varargin{1}*X.value;
                    else
                        error('The first parameter must be of type "double" and the second parameter of type "char"!')
                    end % if
                otherwise
                    error('"unit" can have only one or two parameters')
            end % switch
        end % function
        
        function C=convert(A,B)
            % Method converts a "unit" variable A to the units specified by the string B
            % Syntax:
            %       U=convert(A,'B')
            %                       where A is any variable of class
            %                       "unit" and B is the desired unit. 
            %
            % The value of B may be selected by using the unit.search(str)
            % function or by combining units with operators such as the
            % string 'feet/second'.
            
            if ischar(B)
                % Coerce the type of C to the type of A
                if A.hasBeenConverted
                    A=unit(A);
                end
                C=A;
                % Retain the name of the string in the variable "Name" and
                % create B as a Unit.
                Name=B;
                B=unit(B);
                % Verify that the fundamental units are the same. If they are,
                % make the conversion.
                if (isempty(A.name) && isempty(B.name)) || strcmp(A.name,B.name)
                    % Define the value, offset, and name of C
                    C.value=((A.value+A.offset)-B.offset)/B.value;
                    C.offset=0;
                    C.name=Name;
                    C.hasBeenConverted=true;
                else
                    error(['Incompatible conversion! Fundamental units of ' A.name ' vs. ' B.name])
                end % if
            else
                error('The second parameter in "convert" must be a string')
            end % if
        end % function
        
        % Arithmetic Operators
        function C=uplus(A)
            % The unary "plus" operator for class "unit"
            
            % Unary plus function
            C=A;
        end % function
        
        function C=uminus(A)
            % The unary "minus" operator for class "unit"
            
            % Unary minus function
            C=A;
            C.value=-A.value;
            C.offset=A.offset;
        end % function
        
        function C=plus(A,B)
            % The "plus" function for the addition of class "unit"
            
            % Ensure that A and B are in fundamental units
            A=unit(A);B=unit(B);
            % Check that the units are the same. If not report an error.
%            if A.name==B.name
            if strcmp(A.name,B.name)
                C=A;
                C.value=A.value+B.value;
                C.offset=A.offset+B.offset;
                C.name=A.name;
            else
                error('Incompatible for addition.')
            end % if
        end % function
        
        function C=minus(A,B)
            % The "minus" function for the subtraction of class "unit"
            
            % Ensure that A and B are in fundamental units
            A=unit(A);B=unit(B);
            % Check that the units are the same. If not report an error.
            %if A.name==B.name
            if strcmp(A.name,B.name)
                C=A;
                C.value=A.value-B.value;
                C.offset=A.offset-B.offset;
                C.name=A.name;
            else
                error('Incompatible for subtraction.')
            end % if
        end % function
        
        function C=mtimes(A,B)
            % The "times" function for the matrix multiplication of class "unit"
            
            % If A is a number and B is a "Unit", make C a "Unit" by
            % copying B. Set the value of C from the number A and the value
            % of the "Unit".
            if isa(A,'double') && isa(B,'unit')
                C=B;
                C.value=A*B.value;
                % If A is a "Unit" and b is a number, just permute the
                % operations.
            elseif isa(A,'unit') && isa(B,'double')
                C=A;
                C.value=B*A.value;
            else
                % If neither A or B is a number the type of A and B is not
                % known. Therefore, call "unit" on both to create them as
                % a "Unit". Assign the output variable C the properties of
                % A so that it is known as a "Unit" also.
                A=unit(A);B=unit(B);
                C=A;
                % Assign the value of C as the product of the values of A
                % and B.
                C.value=A.value*B.value;
                % Simplify the units name for both A and B. Retain the
                % structure for each that describes the fundamental (SI)
                % units.
                Ua=simplifyFund(A.name);
                Ub=simplifyFund(B.name);
                % Assign the structure of Ua to a structure for Uc. To
                % multiply the unit strings, add the exponents in the
                % structure. Convert the structure back to a string and
                % assign it to the unit name for C.
                Uc=Ua;
                Uc.exp=Ua.exp+Ub.exp;
                C.name=fundUnits2str(Uc);
            end % if
        end % function
        
        function C=times(A,B)
            % The "times" function for the array multiplication of class "unit"
            
            % If A is a number and B is a "Unit", make C a "Unit" by
            % copying B. Set the value of C from the number A and the value
            % of the "Unit".
            if isa(A,'double') && isa(B,'unit')
                C=B;
                C.value=A.*B.value;
                % If A is a "Unit" and b is a number, just permute the
                % operations.
            elseif isa(A,'unit') && isa(B,'double')
                C=A;
                C.value=B.*A.value;
            else
                % If neither A or B is a number the type of A and B is not
                % known. Therefore, call "unit" on both to create them as
                % a "Unit". Assign the output variable C the properties of
                % A so that it is known as a "Unit" also.
                A=unit(A);B=unit(B);
                C=A;
                % Assign the value of C as the product of the values of A
                % and B.
                C.value=A.value.*B.value;
                % Simplify the units name for both A and B. Retain the
                % structure for each that describes the fundamental (SI)
                % units.
                Ua=simplifyFund(A.name);
                Ub=simplifyFund(B.name);
                % Assign the structure of Ua to a structure for Uc. To
                % multiply the unit strings, add the exponents in the
                % structure. Convert the structure back to a string and
                % assign it to the unit name for C.
                Uc=Ua;
                Uc.exp=Ua.exp+Ub.exp;
                C.name=fundUnits2str(Uc);
            end % if
        end % function
        
        function C=mrdivide(num,div)
            % The "divide" function for the matrix right division of class "unit"
            
            % If the numerator is a "Unit" and the divisor is a number,
            % copy the numerator to C to coerce the type and divide the
            % numerator value by the number.
            if isa(num,'unit') && isa(div,'double')
                C=num;
                C.value=num.value/div;
                % If the numerator is a number and the divisor is a unit,
                % coerce the type of C to the type of the divisor and divide
                % the number in the numerator by the value in the unit.
                % Simplify the name in the divisor and create a structure Udiv.
            elseif isa(div,'unit') && isa(num,'double')
                C=div;
                C.value=num/div.value;
                Udiv=simplifyFund(div.name);
                % Set the output structure equal to the divisor structure.
                Uc=Udiv;
                % Change the sign of the exponent in the divisor to reflect
                % that it is the divisor.
                Uc.exp=-Udiv.exp;
                % Convert the output structure to a string and assign it to
                % the name of the output "Unit".
                C.name=fundUnits2str(Uc);
            else
                % If neither A or B is a number the type of A and B is
                % unknown. Call "unit" on both the numerator and the
                % divisor. Assign the output C to the numerator to coerce
                % the type.
                A=unit(num);B=unit(div);
                C=num;
                % Divide the numerator value by the denominator value and
                % assign it to the output value.
                C.value=A.value/B.value;
                % Simplify both A and B creating a structure for each.
                Ua=simplifyFund(A.name);Ub=simplifyFund(B.name);
                % Create the output structure form A's structure.
                Uc=Ua;
                % Determine the output exponent in the structure by
                % subtracting the divisor exponent from the numerator
                % exponent.
                Uc.exp=Ua.exp-Ub.exp;
                % Convert the output structure to a string and asign it to
                % the output "Unit".
                C.name=fundUnits2str(Uc);
            end % if
        end % function
        
        function C=mldivide(num,div)
            % The "divide" function for matrix left division of class "unit"
            
            % If the numerator is a "Unit" and the divisor is a number,
            % copy the numerator to C to coerce the type and divide the
            % numerator value by the number.
            if isa(num,'unit') && isa(div,'double')
                C=num;
                C.value=num.value\div;
                % If the numerator is a number and the divisor is a unit,
                % coerce the type of C to the type of the divisor and divide
                % the number in the numerator by the value in the unit.
                % Simplify the name in the divisor and create a structure Udiv.
            elseif isa(div,'unit') && isa(num,'double')
                C=div;
                C.value=num\div.value;
                Udiv=simplifyFund(div.name);
                % Set the output structure equal to the divisor structure.
                Uc=Udiv;
                % Change the sign of the exponent in the divisor to reflect
                % that it is the divisor.
                Uc.exp=-Udiv.exp;
                % Convert the output structure to a string and assign it to
                % the name of the output "Unit".
                C.name=fundUnits2str(Uc);
            else
                % If neither A or B is a number the type of A and B is
                % unknown. Call "unit" on both the numerator and the
                % divisor. Assign the output C to the numerator to coerce
                % the type.
                A=unit(num);B=unit(div);
                C=num;
                % Divide the numerator value by the denominator value and
                % assign it to the output value.
                C.value=A.value\B.value;
                % Simplify both A and B creating a structure for each.
                Ua=simplifyFund(A.name);Ub=simplifyFund(B.name);
                % Create the output structure form A's structure.
                Uc=Ua;
                % Determine the output exponent in the structure by
                % subtracting the divisor exponent from the numerator
                % exponent.
                Uc.exp=Ua.exp-Ub.exp;
                % Convert the output structure to a string and asign it to
                % the output "Unit".
                C.name=fundUnits2str(Uc);
            end % if
        end % function
        
        function C=rdivide(num,div)
            % The "divide" function for the arraywise right division of class "unit"
            
            % If the numerator is a "Unit" and the divisor is a number,
            % copy the numerator to C to coerce the type and divide the
            % numerator value by the number.
            if isa(num,'unit') && isa(div,'double')
                C=num;
                C.value=num.value./div;
                % If the numerator is a number and the divisor is a unit,
                % coerce the type of C to the type of the divisor and divide
                % the number in the numerator by the value in the unit.
                % Simplify the name in the divisor and create a structure Udiv.
            elseif isa(div,'unit') && isa(num,'double')
                C=div;
                C.value=num./div.value;
                Udiv=simplifyFund(div.name);
                % Set the output structure equal to the divisor structure.
                Uc=Udiv;
                % Change the sign of the exponent in the divisor to reflect
                % that it is the divisor.
                Uc.exp=-Udiv.exp;
                % Convert the output structure to a string and assign it to
                % the name of the output "Unit".
                C.name=fundUnits2str(Uc);
            else
                % If neither A or B is a number the type of A and B is
                % unknown. Call "unit" on both the numerator and the
                % divisor. Assign the output C to the numerator to coerce
                % the type.
                A=unit(num);B=unit(div);
                C=num;
                % Divide the numerator value by the denominator value and
                % assign it to the output value.
                C.value=A.value./B.value;
                % Simplify both A and B creating a structure for each.
                Ua=simplifyFund(A.name);Ub=simplifyFund(B.name);
                % Create the output structure form A's structure.
                Uc=Ua;
                % Determine the output exponent in the structure by
                % subtracting the divisor exponent from the numerator
                % exponent.
                Uc.exp=Ua.exp-Ub.exp;
                % Convert the output structure to a string and asign it to
                % the output "Unit".
                C.name=fundUnits2str(Uc);
            end % if
        end % function
        
        function C=ldivide(num,div)
            % The "divide" function for the arraywise left division of class "unit"
            
            % If the numerator is a "Unit" and the divisor is a number,
            % copy the numerator to C to coerce the type and divide the
            % numerator value by the number.
            if isa(num,'unit') && isa(div,'double')
                C=num;
                C.value=num.value.\div;
                % If the numerator is a number and the divisor is a unit,
                % coerce the type of C to the type of the divisor and divide
                % the number in the numerator by the value in the unit.
                % Simplify the name in the divisor and create a structure Udiv.
            elseif isa(div,'unit') && isa(num,'double')
                C=div;
                C.value=num.\div.value;
                Udiv=simplifyFund(div.name);
                % Set the output structure equal to the divisor structure.
                Uc=Udiv;
                % Change the sign of the exponent in the divisor to reflect
                % that it is the divisor.
                Uc.exp=-Udiv.exp;
                % Convert the output structure to a string and assign it to
                % the name of the output "Unit".
                C.name=fundUnits2str(Uc);
            else
                % If neither A or B is a number the type of A and B is
                % unknown. Call "Units" on both the numerator and the
                % divisor. Assign the output C to the numerator to coerce
                % the type.
                A=unit(num);B=unit(div);
                C=num;
                % Divide the numerator value by the denominator value and
                % assign it to the output value.
                C.value=A.value.\B.value;
                % Simplify both A and B creating a structure for each.
                Ua=simplifyFund(A.name);Ub=simplifyFund(B.name);
                % Create the output structure form A's structure.
                Uc=Ua;
                % Determine the output exponent in the structure by
                % subtracting the divisor exponent from the numerator
                % exponent.
                Uc.exp=Ua.exp-Ub.exp;
                % Convert the output structure to a string and asign it to
                % the output "Unit".
                C.name=fundUnits2str(Uc);
            end % if
        end % function
        
        function C=mpower(A,B)
            % The "power" function for the matrix power of class "unit"
            
            % Ensure that A is of type "unit" and converted to fundamental
            % units
            A=unit(A);
            % If B is a number, coerce the type of C to a "Unit"
            if isa(B,'double')
                C=A;
                % Set the value of C by raising the value of A to the power
                % of B
                C.value=A.value^B;
                % Simplify A and create a structure
                Ua=simplifyFund(A.name);
                % Raise the fundamental units to the power of B by
                % multiplying the exponents
                Ua.exp=Ua.exp*B;
                % Convert the fundamental units in the structure to a
                % string and assign the string to C.name.
                C.name=fundUnits2str(Ua);
            elseif isa(B,'unit')
                % If B is a "Unit" check that the string is empty so that
                % the operation will make sense.
                if isempty(B.name)
                    % If the string is empty assign A to C to coerce the
                    % type of C
                    C=A;
                    % Use the value of B as an exponent and determine the
                    % value of C.value.
                    exp=B.value;
                    C.value=A.value^exp;
                    % Simplify the name in A and create a structure Ua
                    Ua=simplifyFund(A.name);
                    % Raise the fundamental units to the power of exp
                    Ua.exp=Ua.exp*exp;
                    % Convert the fundamental units in the structure to a
                    % string asn assign the string to C.name.
                    C.name=fundUnits2str(Ua);
                else
                    error('The units string of B must be empty!')
                end % if
            else
                error('An exponent must be a number of type "double"!')
            end % if
        end % function
        
        function C=power(A,B)
            % The "power" function for the array power of class "unit"
            
            % Ensure that A is of type "unit" and converted to fundamental
            % units
            A=unit(A);
            % If B is a number, coerce the type of C to a "Unit"
            if isa(B,'double')
                C=A;
                % Set the value of C by raising the value of A to the power
                % of B
                C.value=A.value.^B;
                % Simplify A and create a structure
                Ua=simplifyFund(A.name);
                % Raise the fundamental units to the power of B by
                % multiplying the exponents
                Ua.exp=Ua.exp*B;
                % Convert the fundamental units in the structure to a
                % string and assign the string to C.name.
                C.name=fundUnits2str(Ua);
            elseif isa(B,'unit')
                % If B is a "Unit" check that the string is empty so that
                % the operation will make sense.
                if isempty(B.name)
                    % If the string is empty assign A to C to coerce the
                    % type of C
                    C=A;
                    % Use the value of B as an exponent and determine the
                    % value of C.value.
                    exp=B.value;
                    C.value=A.value.^exp;
                    % Simplify the name in A and create a structure Ua
                    Ua=simplifyFund(A.name);
                    % Raise the fundamental units to the power of exp
                    Ua.exp=Ua.exp*exp;
                    % Convert the fundamental units in the structure to a
                    % string asn assign the string to C.name.
                    C.name=fundUnits2str(Ua);
                else
                    error('The units string of B must be empty!')
                end % if
            else
                error('An exponent must be a number of type "double"!')
            end % if
        end % function
        
        function C=ctranspose(A)
            % The complex transpose function for class "unit"
            
            % Ensure that A is in fundamental units
            C=unit(A);
            % Transpose the value
            C.value=C.value';
        end %function
        
        function C=transpose(A)
            % The transpose function for class "unit"
            
            % Ensure that A is in fundamental units
            C=unit(A);
            % Transpose the value
            C.value=C.value.';
        end % function
        
        function C=sqrt(A)
            % The square root function for class "unit"
            
            % Ensure that A is in fundamental units
            A=unit(A);
            % Compute the square root of the value
            C=A;
            C.value=sqrt(A.value);
            % Adjust the units
            x=simplifyFund(A.name);
            x.exp=x.exp*0.5;
            C.name=fundUnits2str(x);
        end % function
        
        % Triginometric Functions for Radians
        
        function C=sin(A)
            % The sine function for the class "unit"
            %  The return variable C is type "double"
            
            % Ensure that A is in fundamental units
            A=unit(A);
            % If the name is empty, compute the sin of the value
            if isempty(A.name)
                C=sin(A.value);
            else
                error('The input argument must have no unit name!')
            end % if  
        end % function
        
        function C=csc(A)
            % The cosecant function for the class "unit"
            %  The return variable C is type "double"
            
            % Ensure that A is in fundamental units
            A=unit(A);
            % If the name is empty, compute the csc of the value
            if isempty(A.name)
                C=csc(A.value);
            else
                error('The input argument must have no unit name!')
            end % if  
        end % function
        
        function C=cos(A)
            % The cosine function for the class "unit"
            %  The return variable C is type "double"
            
            % Ensure that A is in fundamental units
            A=unit(A);
            % If the name is empty, compute the cos of the value
            if isempty(A.name)
                C=cos(A.value);
            else
                error('The input argument must have no unit name!')
            end % if  
        end % function
        
        function C=sec(A)
            % The secant function for the class "unit"
            %  The return variable C is type "double"
            
            % Ensure that A is in fundamental units
            A=unit(A);
            % If the name is empty, compute the sec of the value
            if isempty(A.name)
                C=sec(A.value);
            else
                error('The input argument must have no unit name!')
            end % if  
        end % function
        
        function C=tan(A)
            % The tangent function for the class "unit"
            %  The return variable C is type "double"
            
            % Ensure that A is in fundamental units
            A=unit(A);
            % If the name is empty, compute the tan of the value
            if isempty(A.name)
                C=tan(A.value);
            else
                error('The input argument must have no unit name!')
            end % if  
        end % function
        
        function C=cot(A)
            % The cotangent function for the class "unit"
            %  The return variable C is type "double"
            
            % Ensure that A is in fundamental units
            A=unit(A);
            % If the name is empty, compute the cot of the value
            if isempty(A.name)
                C=cot(A.value);
            else
                error('The input argument must have no unit name!')
            end % if  
        end % function
        
        function C=asin(A)
            % The inverse sine function for the class "unit"
            %  The return variable C is type "double"
            
            % Ensure that A is in fundamental units
            A=unit(A);
            % If the name is empty, compute the asin of the value
            if isempty(A.name);
                C=asin(A.value);
            else
                error('The input argument must have no unit name!')
            end
        end % function
        
        function C=acos(A)
            % The inverse cosine function for the class "unit"
            %  The return variable C is type "double"
            
            % Ensure that A is in fundamental units
            A=unit(A);
            % If the name is empty, compute the acos of the value
            if isempty(A.name);
                C=acos(A.value);
            else
                error('The input argument must have no unit name!')
            end
        end % function
        
        function C=atan(A)
            % The inverse tangent function for the class "unit"
            %  The return variable C is type "double"
            
            % Ensure that A is in fundamental units
            A=unit(A);
            % If the name is empty, compute the atan of the value
            if isempty(A.name);
                C=atan(A.value);
            else
                error('The input argument must have no unit name!')
            end
        end % function
        
        function C=acot(A)
            % The inverse cotangent function for the class "unit"
            %  The return variable C is type "double"
            
            % Ensure that A is in fundamental units
            A=unit(A);
            % If the name is empty, compute the acot of the value
            if isempty(A.name);
                C=acot(A.value);
            else
                error('The input argument must have no unit name!')
            end
        end % function
        
        function C=asec(A)
            % The inverse secant function for the class "unit"
            %  The return variable C is type "double"
            
            % Ensure that A is in fundamental units
            A=unit(A);
            % If the name is empty, compute the asec of the value
            if isempty(A.name);
                C=asec(A.value);
            else
                error('The input argument must have no unit name!')
            end
        end % function
        
        function C=acsc(A)
            % The inverse cosecant function for the class "unit"
            %  The return variable C is type "double"
            
            % Ensure that A is in fundamental units
            A=unit(A);
            % If the name is empty, compute the acsc of the value
            if isempty(A.name);
                C=acsc(A.value);
            else
                error('The input argument must have no unit name!')
            end
        end % function
        
        % Hyperbolic functions for radians
        
        function C=sinh(A)
            % The hyperbolic sine function for the class "unit"
            %  The return variable C is type "double"
            
            % Ensure that A is in fundamental units
            A=unit(A);
            % If the name is empty, compute the sinh of the value
            if isempty(A.name)
                C=sinh(A.value);
            else
                error('The input argument must have no unit name!')
            end % if  
        end % function
        
        function C=cosh(A)
            % The hyperbolic cosine function for the class "unit"
            %  The return variable C is type "double"
            
            % Ensure that A is in fundamental units
            A=unit(A);
            % If the name is empty, compute the cosh of the value
            if isempty(A.name)
                C=cosh(A.value);
            else
                error('The input argument must have no unit name!')
            end % if  
        end % function
        
        function C=tanh(A)
            % The hyperbolic tangent function for the class "unit"
            %  The return variable C is type "double"
            
            % Ensure that A is in fundamental units
            A=unit(A);
            % If the name is empty, compute the tanh of the value
            if isempty(A.name)
                C=tanh(A.value);
            else
                error('The input argument must have no unit name!')
            end % if  
        end % function
        
        function C=csch(A)
            % The hyperbolic cosecant function for the class "unit"
            %  The return variable C is type "double"
            
            % Ensure that A is in fundamental units
            A=unit(A);
            % If the name is empty, compute the csch of the value
            if isempty(A.name)
                C=csch(A.value);
            else
                error('The input argument must have no unit name!')
            end % if  
        end % function
        
        function C=sech(A)
            % The hyperbolic secant function for the class "unit"
            %  The return variable C is type "double"
            
            % Ensure that A is in fundamental units
            A=unit(A);
            % If the name is empty, compute the sech of the value
            if isempty(A.name)
                C=sech(A.value);
            else
                error('The input argument must have no unit name!')
            end % if  
        end % function
        
        function C=coth(A)
            % The hyperbolic cotangent function for the class "unit"
            %  The return variable C is type "double"
            
            % Ensure that A is in fundamental units
            A=unit(A);
            % If the name is empty, compute the coth of the value
            if isempty(A.name)
                C=coth(A.value);
            else
                error('The input argument must have no unit name!')
            end % if  
        end % function
        
        function C=asinh(A)
            % The inverse hyperbolic sine function for the class "unit"
            %  The return variable C is type "double"
            
            % Ensure that A is in fundamental units
            A=unit(A);
            % If the name is empty, compute the asinh of the value
            if isempty(A.name);
                C=asinh(A.value);
            else
                error('The input argument must have no unit name!')
            end
        end % function
        
        function C=acosh(A)
            % The inverse hyperbolic cosine function for the class "unit"
            %  The return variable C is type "double"
            
            % Ensure that A is in fundamental units
            A=unit(A);
            % If the name is empty, compute the acosh of the value
            if isempty(A.name);
                C=acosh(A.value);
            else
                error('The input argument must have no unit name!')
            end
        end % function
        
        function C=atanh(A)
            % The inverse hyperbolic tangent for the class "unit"
            %  The return variable C is type "double"
            
            % Ensure that A is in fundamental units
            A=unit(A);
            % If the name is empty, compute the atanh of the value
            if isempty(A.name);
                C=atanh(A.value);
            else
                error('The input argument must have no unit name!')
            end
        end % function
        
        function C=acoth(A)
            % The inverse hyperbolic cotangent for the class "unit"
            
            % Ensure that A is in fundamental units
            A=unit(A);
            % If the name is empty, compute the acoth of the value
            if isempty(A.name);
                C=acoth(A.value);
            else
                error('The input argument must have no unit name!')
            end
        end % function
        
        function C=asech(A)
            % The inverse hyperbolic secant function for the class "unit"
            %  The return variable C is type "double"
            
            % Ensure that A is in fundamental units
            A=unit(A);
            % If the name is empty, compute the asech of the value
            if isempty(A.name);
                C=asech(A.value);
            else
                error('The input argument must have no unit name!')
            end
        end % function
        
        function C=acsch(A)
            % The inverse hyperbolic cosecant function for the class "unit"
            %  The return variable C is type "double"
            
            % Ensure that A is in fundamental units
            A=unit(A);
            % If the name is empty, compute the acsch of the value
            if isempty(A.name);
                C=acsch(A.value);
            else
                error('The input argument must have no unit name!')
            end
        end % function
        
        % Exponential Functions
        function C=log10(A)
            % The function for the common logrithm of the class "unit
            %  The return variable C is type "double"
            
            % Ensure that A has fundamental units
            A=unit(A);
            % Take the log of the value and return it in C
            C=log10(A.value+A.offset);
        end % function
        
        function C=log(A)
            % The function for the natural logrithm of the class "unit"
            %  The return variable C is type "double"
            
            % Ensure that A has fundamental units
            A=unit(A);
            % Take the natural log of the value and return it in C
            C=log(A.value+A.offset);
        end % function
        
        function C=exp(A)
            % The function raising "e" to a power of the class "unit"
            %  The return variable C is type "double"
            
            % Ensure that A is in fundamental units
            A=unit(A);
            % Raise the value to the power of "e"
            C=exp(A.value);
        end % function
        
        % Logical Functions
        function C=eq(A,B)
            % Compares unit A and B for equaltiy
            
            % It is not known what the class of A and B is, so call "unit"
            % on both A and B.
            Ua=unit(A);
            Ub=unit(B);
            % If the offset of Ua is nonzero, add value and offset before
            % comparison
            if Ua.offset~=0
                Ua.value=Ua.value+Ua.offset;
                Ua.offset=0;
            end % if
            % If the offset of Ub is nonzero, add value and offset before
            % comparison
            if Ub.offset~=0
                Ub.value=Ub.value+Ub.offset;
                Ub.offset=0;
            end % if
            % Compare the two values for A and B by subtracting them and
            % taking the absolute value of the difference. Use "eps" to
            % find the next larger floating point number for the Ua and Ub
            % value. Double the sum of these values and call the values
            % equal if the absolute value of the difference is less than
            % double the sum.
            L1=abs(Ua.value-Ub.value)<2*(eps(Ua.value)+eps(Ub.value));
            % Simplify the name string for Ua and Ub and create a structure
            % for each.
            Fa=simplifyFund(Ua.name);
            Fb=simplifyFund(Ub.name);
            % Test that the exponents in the name strin are equal
            L2=eq(Fa.exp,Fb.exp);
            % Test that the values are equal and that all the exponents in
            % the name strings are equal. If true set C to be true.
            % Otherwise, set C to be false.
            if all(all(L1)) && all(L2)
                C=true;
            else
                C=false;
            end % if
        end % function
        
        function C=ne(A,B)
            % Compares unit A nd B to determine if they are not equal
            
            % It is not known what the class of A and B is, so call "unit"
            % on both A and B.
            Ua=unit(A);
            Ub=unit(B);
            % If the offset of Ua is nonzero, add value and offset before
            % comparison
            if Ua.offset~=0
                Ua.value=Ua.value+Ua.offset;
                Ua.offset=0;
            end % if
            % If the offset of Ub is nonzero, add value and offset before
            % comparison
            if Ub.offset~=0
                Ub.value=Ub.value+Ub.offset;
                Ub.offset=0;
            end % if
            % Compare the two values for A and B by subtracting them and
            % taking the absolute value of the difference. Use "eps" to
            % find the next larger floating point number for the Ua and Ub
            % value. Double the sum of these values and call the values
            % equal if the absolute value of the difference is less than
            % double the sum.
            L1=abs(Ua.value-Ub.value)>=2*(eps(Ua.value)+eps(Ub.value));
            % Simplify the name string for Ua and Ub and create a structure
            % for each.
            Fa=simplifyFund(Ua.name);
            Fb=simplifyFund(Ub.name);
            % Test that the exponents in the name string are not equal
            L2=ne(Fa.exp,Fb.exp);
            % Test that any of the values are not equal and that all the
            % exponents in the name strings are equal. If true, set C to be
            % true. Otherwise, set C to be false.
            if all(all(L1)) || any(L2)
                C=true;
            else
                C=false;
            end % if
        end % function
        
        function C=lt(A,B)
            % Compares unit A to B to determine if A is less than B
            
            % It is not known what the class of A and B is, so call "unit"
            % on both A and B.
            Ua=unit(A);
            Ub=unit(B);
            % If the offset of Ua is nonzero, add value and offset before
            % comparison
            if Ua.offset~=0
                Ua.value=Ua.value+Ua.offset;
                Ua.offset=0;
            end % if
            % If the offset of Ub is nonzero, add value and offset before
            % comparison
            if Ub.offset~=0
                Ub.value=Ub.value+Ub.offset;
                Ub.offset=0;
            end % if
            % Compare the values
            L1=Ua.value<Ub.value;
            % Simplify the name string for Ua and Ub and create a structure
            % for each.
            Fa=simplifyFund(Ua.name);
            Fb=simplifyFund(Ub.name);
            % Test that the exponents in the name string are equal
            L2=eq(Fa.exp,Fb.exp);
            % Test that any of the values are not equal and that all the
            % exponents in the name strings are equal. If true, set C to be
            % true. Otherwise, set C to be false.
            if all(all(L1)) && all(L2)
                C=true;
            else
                C=false;
            end % if
        end % function
        
        function C=le(A,B)
            % Compares unit to determine if A is less than or equal to B
            
            % It is not known what the class of A and B is, so call "unit"
            % on both A and B.
            Ua=unit(A);
            Ub=unit(B);
            % If the offset of Ua is nonzero, add value and offset before
            % comparison
            if Ua.offset~=0
                Ua.value=Ua.value+Ua.offset;
                Ua.offset=0;
            end % if
            % If the offset of Ub is nonzero, add value and offset before
            % comparison
            if Ub.offset~=0
                Ub.value=Ub.value+Ub.offset;
                Ub.offset=0;
            end % if
            % Compare the values
            L1=Ua.value<=Ub.value;
            % Simplify the name string for Ua and Ub and create a structure
            % for each.
            Fa=simplifyFund(Ua.name);
            Fb=simplifyFund(Ub.name);
            % Test that the exponents in the name string are equal
            L2=eq(Fa.exp,Fb.exp);
            % Test that any of the values are not equal and that all the
            % exponents in the name strings are equal. If true, set C to be
            % true. Otherwise, set C to be false.
            if all(all(L1)) && all(L2)
                C=true;
            else
                C=false;
            end % if
        end % function
        
        function C=gt(A,B)
            % Compares unit A to B to determine if A is greater than B
            
            % It is not known what the class of A and B is, so call "unit"
            % on both A and B.
            Ua=unit(A);
            Ub=unit(B);
            % If the offset of Ua is nonzero, add value and offset before
            % comparison
            if Ua.offset~=0
                Ua.value=Ua.value+Ua.offset;
                Ua.offset=0;
            end % if
            % If the offset of Ub is nonzero, add value and offset before
            % comparison
            if Ub.offset~=0
                Ub.value=Ub.value+Ub.offset;
                Ub.offset=0;
            end % if
            % Compare the values
            L1=Ua.value>Ub.value;
            % Simplify the name string for Ua and Ub and create a structure
            % for each.
            Fa=simplifyFund(Ua.name);
            Fb=simplifyFund(Ub.name);
            % Test that the exponents in the name string are equal
            L2=eq(Fa.exp,Fb.exp);
            % Test that any of the values are not equal and that all the
            % exponents in the name strings are equal. If true, set C to be
            % true. Otherwise, set C to be false.
            if all(all(L1)) && all(L2)
                C=true;
            else
                C=false;
            end % if
        end % Function
        
        function C=ge(A,B)
            % Compares unit to determine if A is greater than or equal to B
            
            % It is not known what the class of A and B is, so call "unit"
            % on both A and B.
            Ua=unit(A);
            Ub=unit(B);
            % If the offset of Ua is nonzero, add value and offset before
            % comparison
            if Ua.offset~=0
                Ua.value=Ua.value+Ua.offset;
                Ua.offset=0;
            end % if
            % If the offset of Ub is nonzero, add value and offset before
            % comparison
            if Ub.offset~=0
                Ub.value=Ub.value+Ub.offset;
                Ub.offset=0;
            end % if
            % Compare the values
            L1=Ua.value>=Ub.value;
            % Simplify the name string for Ua and Ub and create a structure
            % for each.
            Fa=simplifyFund(Ua.name);
            Fb=simplifyFund(Ub.name);
            % Test that the exponents in the name string are equal
            L2=eq(Fa.exp,Fb.exp);
            % Test that any of the values are not equal and that all the
            % exponents in the name strings are equal. If true, set C to be
            % true. Otherwise, set C to be false.
            if all(all(L1)) && all(L2)
                C=true;
            else
                C=false;
            end % if
        end % function
        
    % Plotting methods
        function [handles,legendHandle]=plot(varargin)
        % This method plots an X and Y vector where X and Y are "unit"
        % SYNTAX:
        %   plot(X,Y)
        %       where X must be a "unit" vector and Y can be either a
        %       "unit" vector or a structure of "unit" vectors. The
        %       vector length of Y must be the same as the vector length of
        %       X. The default labels for the plot are the same as the
        %       variable names used in the call to "unit". Also included
        %       in the axis labels are the units of the X and Y variables.
        %       If the Y variable is a structure where each field contains
        %       a "unit" vector (the same length as X), then all fields
        %       are plotted and the field names are used to create the
        %       legend.
        %   plot(X,Y,LineSpec) where LineSpec is the same as the standard
        %       plot statement.
        %   plot(X,Y,LineSpec,'PropertyName',PropertyValue,...)
        %       Where PropertyName and PropertyValue are the same pairs
        %       used in the standard plot statement.
        %   h=plot(X,Y,LineSpec,'PropertyName',PropertyValue,...) returns
        %       the handles to the lineseries objects, one handle per line.
        
        % Find the names of the input variables
        count=1;
        for I=1:nargin
            temp=inputname(I);
            if ~isempty(temp)
                varname{count,1}=temp;
            end %if
            count=count+1;
        end % for
        
        % Create figure
        figure1 = figure;
        % Create axes
        axes1 = axes('Parent',figure1);
        box(axes1,'on');
        
        % BEGIN PLOTTING       
        % Plot the case where there are two "unit" vectors
        [rows,~]=size(varname);       
        if rows==2 && isa(varargin{1},'unit') && isa(varargin{2},'unit')
            X=varargin{1};
            Y=varargin{2};
            % If there are more arguments than two use them
            if nargin>2
                h=plot(X.value+X.offset,Y.value+Y.offset,varargin{3:nargin});
            else
                h=plot(X.value+X.offset,Y.value+Y.offset);
            end
            % Create the axis labels and use them
            xAxisLabel=[varname{1} ' (' X.name ')'];
            yAxisLabel=[varname{2} ' (' Y.name ')'];
            xlabel(xAxisLabel);
            ylabel(yAxisLabel);
            
        % Plot the case where there is a "unit" vector for X and a unit structure
        % for Y. The fields in Y must each have the same vector length as the as
        % the fields in X.
        elseif rows==2 && isa(varargin{1},'unit') && isa(varargin{2},'struct')
            X=varargin{1};
            temp=varargin{2};
            Y=[];
            % Obtain the field names
            lineNames=fieldnames(varargin{2});
            L=length(X.value);
            N=numel(lineNames);
            % Check that the new row in the matrix of Y all have the same
            % length as the X variable and then create the Y matrix
            for I=1:N
                if isa(temp.(lineNames{I}),'unit')
                    newRow=temp.(lineNames{I}).value+temp.(lineNames{I}).offset;
                    if length(newRow)~=L;
                        error('All rows in a "unit" must be the same length as the X value')
                    end % if
                    [row,col]=size(newRow);
                    if col>row
                        Y=[Y;newRow];
                    else
                        Y=[Y;newRow.'];
                    end %if
                else
                    error(['All fields of ' varname(2) ' must be of type "unit"'])
                end % if
            end % for
            % If there are more arguments in the plot statement, use them
            if nargin>2
                h=plot(X.value+X.offset,Y,varargin{3:nargin});
            else
                h=plot(X.value+X.offset,Y);
            end
            % Create the axis labels and use them
            xAxisLabel=[varname{1} ' (' X.name ')'];
            yAxisLabel=[varname{2} ' (' temp.(lineNames{1}).name ')'];
            xlabel(xAxisLabel);
            ylabel(yAxisLabel);
            % Create the legend names from the field names
            for I=1:N
                set(h(I),'DisplayName',lineNames{I})
            end % for
            % Enable the legend
            lh=legend(axes1,'show');
        else
            error('Cannot not plot the "y" variable perhaps because it is a structure and field. Change to a simple variable name.')
        end % if
        % If an output argument was requested, construct it
        if nargout>0
            handles=h;
            if isa(varargin{2},'struct')
                legendHandle=lh;
            end %if
        end
    end % function
    
        function [handles,legendHandle]=semilogx(varargin)
        % This method plots an X and Y vector where X and Y are "unit"
        % SYNTAX:
        %   semilogx(X,Y)
        %       where X must be a "unit" vector and Y can be either a
        %       "unit" vector or a structure of "unit" vectors. The
        %       vector length of Y must be the same as the vector length of
        %       X. The default labels for the semilogx are the same as the
        %       variable names used in the call to "unit". Also included
        %       in the axis labels are the units of the X and Y variables.
        %       If the Y variable is a structure where each field contains
        %       a "unit" vector (the same length as X), then all fields
        %       are plotted and the field names are used to create the
        %       legend.
        %   semilogx(X,Y,LineSpec) where LineSpec is the same as the standard
        %       semilogx statement.
        %   semilogx(X,Y,LineSpec,'PropertyName',PropertyValue,...)
        %       Where PropertyName and PropertyValue are the same pairs
        %       used in the standard semilogx statement.
        %   h=semilogx(X,Y,LineSpec,'PropertyName',PropertyValue,...) returns
        %       the handles to the lineseries objects, one handle per line.
        
        % Find the names of the input variables
        count=1;
        for I=1:nargin
            temp=inputname(I);
            if ~isempty(temp)
                varname{count,1}=temp;
            end %if
            count=count+1;
        end % for
        
        % Create figure
        figure1 = figure;
        % Create axes
        axes1 = axes('Parent',figure1);
        box(axes1,'on');
        
        % BEGIN PLOTTING       
        % Plot the case where there are two "unit" vectors
        [rows,~]=size(varname);  
        if rows==2 && isa(varargin{1},'unit') && isa(varargin{2},'unit')
            X=varargin{1};
            Y=varargin{2};
            % If there are more arguments than two use them
            if nargin>2
                h=semilogx(X.value+X.offset,Y.value+Y.offset,varargin{3:nargin});
            else
                h=semilogx(X.value+X.offset,Y.value+Y.offset);
            end
            % Create the axis labels and use them
            xAxisLabel=[varname{1} ' (' X.name ')'];
            yAxisLabel=[varname{2} ' (' Y.name ')'];
            xlabel(xAxisLabel);
            ylabel(yAxisLabel);
            
        % Plot the case where there is a "unit" vector for X and a unit structure
        % for Y. The fields in Y must each have the same vector length as the as
        % the fields in X.
        elseif rows==2 && isa(varargin{1},'unit') && isa(varargin{2},'struct')
            X=varargin{1};
            temp=varargin{2};
            Y=[];
            % Obtain the field names
            lineNames=fieldnames(varargin{2});
            L=length(X.value);
            N=numel(lineNames);
            % Check that the new row in the matrix of Y all have the same
            % length as the X variable and then create the Y matrix
            for I=1:N
                if isa(temp.(lineNames{I}),'unit')
                    newRow=temp.(lineNames{I}).value+temp.(lineNames{I}).offset;
                    if length(newRow)~=L;
                        error('All rows in a "unit" must be the same length as the X value')
                    end % if
                    [row,col]=size(newRow);
                    if col>row
                        Y=[Y;newRow];
                    else
                        Y=[Y;newRow.'];
                    end %if
                else
                    error(['All fields of ' varname(2) ' must be of type "unit"'])
                end % if
            end % for
            % If there are more arguments in the semilogx statement, use them
            if nargin>2
                h=semilogx(X.value+X.offset,Y,varargin{3:nargin});
            else
                h=semilogx(X.value+X.offset,Y);
            end
            % Create the axis labels and use them
            xAxisLabel=[varname{1} ' (' X.name ')'];
            yAxisLabel=[varname{2} ' (' temp.(lineNames{1}).name ')'];
            xlabel(xAxisLabel);
            ylabel(yAxisLabel);
            % Create the legend names from the field names
            for I=1:N
                set(h(I),'DisplayName',lineNames{I})
            end % for
            % Enable the legend
            lh=legend(axes1,'show');
        else
            error('Cannot not plot the "y" variable perhaps because it is a structure and field. Change to a simple variable name.')
        end % if
        % If an output argument was requested, construct it
        if nargout>0
            handles=h;
            legendHandle=lh;
        end
    end % function
        
        function [handles,legendHandle]=semilogy(varargin)
        % This method plots an X and Y vector where X and Y are "unit"
        % SYNTAX:
        %   semilogy(X,Y)
        %       where X must be a "unit" vector and Y can be either a
        %       "unit" vector or a structure of "unit" vectors. The
        %       vector length of Y must be the same as the vector length of
        %       X. The default labels for the semilogy are the same as the
        %       variable names used in the call to "unit". Also included
        %       in the axis labels are the units of the X and Y variables.
        %       If the Y variable is a structure where each field contains
        %       a "unit" vector (the same length as X), then all fields
        %       are plotted and the field names are used to create the
        %       legend.
        %   semilogy(X,Y,LineSpec) where LineSpec is the same as the standard
        %       semilogy statement.
        %   semilogy(X,Y,LineSpec,'PropertyName',PropertyValue,...)
        %       Where PropertyName and PropertyValue are the same pairs
        %       used in the standard semilogy statement.
        %   h=semilogy(X,Y,LineSpec,'PropertyName',PropertyValue,...) returns
        %       the handles to the lineseries objects, one handle per line.
        
        % Find the names of the input variables
        count=1;
        for I=1:nargin
            temp=inputname(I);
            if ~isempty(temp)
                varname{count,1}=temp;
            end %if
            count=count+1;
        end % for
        
        % Create figure
        figure1 = figure;
        % Create axes
        axes1 = axes('Parent',figure1);
        box(axes1,'on');
        
        % BEGIN PLOTTING       
        % Plot the case where there are two "unit" vectors
        [rows,~]=size(varname); 
        if rows==2 && isa(varargin{1},'unit') && isa(varargin{2},'unit')
            X=varargin{1};
            Y=varargin{2};
            % If there are more arguments than two use them
            if nargin>2
                h=semilogy(X.value+X.offset,Y.value+Y.offset,varargin{3:nargin});
            else
                h=semilogy(X.value+X.offset,Y.value+Y.offset);
            end
            % Create the axis labels and use them
            xAxisLabel=[varname{1} ' (' X.name ')'];
            yAxisLabel=[varname{2} ' (' Y.name ')'];
            xlabel(xAxisLabel);
            ylabel(yAxisLabel);
            
        % Plot the case where there is a "unit" vector for X and a unit structure
        % for Y. The fields in Y must each have the same vector length as the as
        % the fields in X.
        elseif rows==2 && isa(varargin{1},'unit') && isa(varargin{2},'struct')
            X=varargin{1};
            temp=varargin{2};
            Y=[];
            % Obtain the field names
            lineNames=fieldnames(varargin{2});
            L=length(X.value);
            N=numel(lineNames);
            % Check that the new row in the matrix of Y all have the same
            % length as the X variable and then create the Y matrix
            for I=1:N
                if isa(temp.(lineNames{I}),'unit')
                    newRow=temp.(lineNames{I}).value+temp.(lineNames{I}).offset;
                    if length(newRow)~=L;
                        error('All rows in a "unit" must be the same length as the X value')
                    end % if
                    [row,col]=size(newRow);
                    if col>row
                        Y=[Y;newRow];
                    else
                        Y=[Y;newRow.'];
                    end %if
                else
                    error(['All fields of ' varname(2) ' must be of type "unit"'])
                end % if
            end % for
            % If there are more arguments in the semilogy statement, use them
            if nargin>2
                h=semilogy(X.value+X.offset,Y,varargin{3:nargin});
            else
                h=semilogy(X.value+X.offset,Y);
            end
            % Create the axis labels and use them
            xAxisLabel=[varname{1} ' (' X.name ')'];
            yAxisLabel=[varname{2} ' (' temp.(lineNames{1}).name ')'];
            xlabel(xAxisLabel);
            ylabel(yAxisLabel);
            % Create the legend names from the field names
            for I=1:N
                set(h(I),'DisplayName',lineNames{I})
            end % for
            % Enable the legend
            lh=legend(axes1,'show');
        else
            error('Cannot not plot the "y" variable perhaps because it is a structure and field. Change to a simple variable name.')
        end % if
        % If an output argument was requested, construct it
        if nargout>0
            handles=h;
            legendHandle=lh;
        end
    end % function
    
        function [handles,legendHandle]=loglog(varargin)
        % This method plots an X and Y vector where X and Y are "unit"
        % SYNTAX:
        %   loglog(X,Y)
        %       where X must be a "unit" vector and Y can be either a
        %       "unit" vector or a structure of "unit" vectors. The
        %       vector length of Y must be the same as the vector length of
        %       X. The default labels for the loglog are the same as the
        %       variable names used in the call to "unit". Also included
        %       in the axis labels are the units of the X and Y variables.
        %       If the Y variable is a structure where each field contains
        %       a "unit" vector (the same length as X), then all fields
        %       are plotted and the field names are used to create the
        %       legend.
        %   loglog(X,Y,LineSpec) where LineSpec is the same as the standard
        %       loglog statement.
        %   loglog(X,Y,LineSpec,'PropertyName',PropertyValue,...)
        %       Where PropertyName and PropertyValue are the same pairs
        %       used in the standard loglog statement.
        %   h=loglog(X,Y,LineSpec,'PropertyName',PropertyValue,...) returns
        %       the handles to the lineseries objects, one handle per line.
        
        % Find the names of the input variables
        count=1;
        for I=1:nargin
            temp=inputname(I);
            if ~isempty(temp)
                varname{count,1}=temp;
            end %if
            count=count+1;
        end % for
        
        % Create figure
        figure1 = figure;
        % Create axes
        axes1 = axes('Parent',figure1);
        box(axes1,'on');
        
        % BEGIN PLOTTING       
        % Plot the case where there are two "unit" vectors
        [rows,~]=size(varname); 
        if rows==2 && isa(varargin{1},'unit') && isa(varargin{2},'unit')
            X=varargin{1};
            Y=varargin{2};
            % If there are more arguments than two use them
            if nargin>2
                h=loglog(X.value+X.offset,Y.value+Y.offset,varargin{3:nargin});
            else
                h=loglog(X.value+X.offset,Y.value+Y.offset);
            end
            % Create the axis labels and use them
            xAxisLabel=[varname{1} ' (' X.name ')'];
            yAxisLabel=[varname{2} ' (' Y.name ')'];
            xlabel(xAxisLabel);
            ylabel(yAxisLabel);
            
        % Plot the case where there is a "unit" vector for X and a unit structure
        % for Y. The fields in Y must each have the same vector length as the as
        % the fields in X.
        elseif rows==2 && isa(varargin{1},'unit') && isa(varargin{2},'struct')
            X=varargin{1};
            temp=varargin{2};
            Y=[];
            % Obtain the field names
            lineNames=fieldnames(varargin{2});
            L=length(X.value);
            N=numel(lineNames);
            % Check that the new row in the matrix of Y all have the same
            % length as the X variable and then create the Y matrix
            for I=1:N
                if isa(temp.(lineNames{I}),'unit')
                    newRow=temp.(lineNames{I}).value+temp.(lineNames{I}).offset;
                    if length(newRow)~=L;
                        error('All rows in a "unit" must be the same length as the X value')
                    end % if
                    [row,col]=size(newRow);
                    if col>row
                        Y=[Y;newRow];
                    else
                        Y=[Y;newRow.'];
                    end %if
                else
                    error(['All fields of ' varname(2) ' must be of type "unit"'])
                end % if
            end % for
            % If there are more arguments in the loglog statement, use them
            if nargin>2
                h=loglog(X.value+X.offset,Y,varargin{3:nargin});
            else
                h=loglog(X.value+X.offset,Y);
            end
            % Create the axis labels and use them
            xAxisLabel=[varname{1} ' (' X.name ')'];
            yAxisLabel=[varname{2} ' (' temp.(lineNames{1}).name ')'];
            xlabel(xAxisLabel);
            ylabel(yAxisLabel);
            % Create the legend names from the field names
            for I=1:N
                set(h(I),'DisplayName',lineNames{I})
            end % for
            % Enable the legend
            lh=legend(axes1,'show');
        else
            error('Cannot not plot the "y" variable perhaps because it is a structure and field. Change to a simple variable name.')
        end % if
        % If an output argument was requested, construct it
        if nargout>0
            handles=h;
            legendHandle=lh;
        end
    end % function  
        
    % Utility Functions    
        function display(x)
            % This method displays a unit class on the screen
            %disp(inputname(1))
            if numel(x.value)==1
                disp(' ')
                disp([inputname(1) ' = '])
                s=sprintf('     %g %s',x.value+x.offset,x.name);
                disp(s)
            else
                disp(' ')
                disp([inputname(1) ' = '])
                s=sprintf('    unit of  "%s"',x.name);
                disp(s)
                disp(x.value+x.offset)
            end %if
        end % function
                
    end % methods
    
    methods (Static)
        
        function search(string)
            % Find all unit abbreviations that pertain to the input string
            % Syntax:
            %       unit.search(string)
            %
            % This method searches through the comments relating to all units
            % to find those abbreviations that contain the string
            
            % Convert all strings to lower case
            string=lower(string);
            % Searches for the "string" in the case comment field
            [abbr,comments]=unit.available_units;
            N=numel(abbr);
            count=0;
            % Search the comment strings. If the string is "all" return all
            % unit names
            if strcmp(string,'all')
                index=1:numel(comments);
            else
                for I=1:N
                    stringExists=strfind(comments{I},string);
                    if ~isempty(stringExists{1})
                        count=count+1;
                        index(count)=I;
                    end %if
                end % for
            end %if
            % Display the results
            for I=1:numel(index)
                s=sprintf('     Units of %s     %s',char(comments{index(I)}),char(abbr{index(I)}));
                disp(s)
            end; % for
        end % function
        
        function [unitType,label]=available_units
            % Reads the "define" method to get unit abreviations
            % Syntax:
            %       [unitType,label]=unit.available_units
            %           where unitType and label are cells
            
            % Open the file containing the "define" method
            fid=fopen('convConstant.m');
            % Read each line of the file into a cell array
            file=textscan(fid,'%s','Delimiter','\n');
            % Close the file
            fclose(fid);
            % Simplify the cell array into a simple cell list
            list=file{:};
            % Find each line that begins with "case" and record the index in "list"
            % where it is found. Lines without a beginning "case" will have an index
            % of zero.
            count=0;
            N=numel(list);
            caseIndex(N)=0;
            for I=1:N
                count=count+1;
                line=strfind(list{I},'case');
                if ~isempty(line)
                    caseIndex(count)=I;
                end % if
            end % for
            % Remove all zero values in the "caseIndex" and retain only the lines from
            % file that begin with "case"
            caseIndex= caseIndex>0;
            caseLines=file{:}(caseIndex);
            % Split the lines at the "%" symbol. Trailing comments are in the second
            % part of "A"
            N=numel(caseLines);
            A{N}=[];
            for I=1:N
                A{I}=textscan(caseLines{I},'%s %s','Delimiter','%');
            end % for
            % Place each comment in in a cell list called "label"
            label{N}=[];
            for I=1:N
                label{I}=lower(A{I}{2});
            end % for
            % From the first part of "A" remove the "case" characters and then trim all
            % leading and trailing blanks from the string. The result is the list of
            % permissable units
            unitType{N}=[];
            for I=1:N
                B=A{I}{1};
                B=strrep(B,'case','');
                B=strtrim(B);
                unitType{I}=B;
            end % for
        end %function
        
    end % Static methods
    
end % classdef

function fundUnit=simplifyFund(s)
% Converts a string of fundamental units to a structure.
% Converts a string containing only fundamental units combined with "1",
% "*", "/" and "^" to a structure with names and associated exponents.
% For example, if
%   s='kg*kg*K^-4/m^2'
% returns a structure containing
%   fundUnit =
%              name: {'m'  'kg'  's'  'A'  'K'  'mol'  'cd'}
%              exp: [-2 2 0 0 -4 0 0]
% USEAGE:
%   fundUnit=simplifyFund(s)

% Define a structure to hold fundamental units and there associated
% exponents
fundUnit.name={'m','kg','s','A','K','mol','cd'};
f=length(fundUnit.name);
fundUnit.exp=zeros(1,f);

% Initialize the position index to the first letter of the string "s"
p=1;

% Determine the length of the input string
L=length(s);

if L>0
    % Initialize the value of the token that will hold a fundamental unit
    token='';
    
    % Initialize the value of the increment that changes the exponent. The
    % permissable values of the increment are 1, 0, or -1
    inc=1;
    
    % Test to see if the first letter of the string is alphanumeric
    if ~isstrprop(s(p),'alphanum')
        error('The first character of a unit string must be "1" of a letter');
    end
    % Examine each letter in the string
    while p<=L
        % If this character is a letter, place it in the token and move to the
        % next character
        if isletter(s(p))
            token=[token s(p)]; 
            p=p+1;
            % If the character is not a letter but it is the first character check
            % to see if it is a 1 and move the index to the next character
        elseif isstrprop(s(p),'digit')&&p==1
            if s(p)=='1'
                p=p+1;
            else
                error('A leading numeric character in a units string must be a "1"!')
            end
            % If the character is not a letter or the number "1", find the index of
            % the fundamental unit represented by the token
        else
            ind=strcmp(token,fundUnit.name);
            % Check to be sure that the token is a fundamental unit
            if ~any(ind)==1 && ~isempty(token)
                error(['The string "' token '" is not a fundamental unit!'])
            end %if
            % Find the operator "*","/", or "^"
            switch s(p)
                case '*'
                    % For the token to the left of the "*", add the increment
                    % to the exponent, set the increment for the comming token
                    % (on the right) and clear the current token
                    fundUnit.exp(ind)=fundUnit.exp(ind)+inc;
                    inc=1;
                    token='';
                case '/'
                    % For the token to the leftof the "/", add the increment to
                    % the exponent, set the increment for the comming token on
                    % the right to "-1" and clear the current token
                    fundUnit.exp(ind)=fundUnit.exp(ind)+inc;
                    inc=-1;
                    token='';
                case '^'
                    % Create a string to hold the exponent
                    digstr='';
                    % If the current character is a "-" sign and we are not at
                    % the end of the input string, place the minus sign in the
                    % digit string and move to the next character
                    if p<L && s(p+1)=='-'
                        digstr=[digstr s(p+1)]; 
                        p=p+1;
                    end %if
                    % If we at least two characters the end and the next two
                    % characters are digits, add them to the digit string and
                    % advance the position pointer by two characters
                    if p+2<=L && all(isstrprop(s(p+1:p+2),'digit'))
                        digstr=[digstr s(p+1:p+2)]; 
                        p=p+2;
                        % If we are at least one character from the end and the
                        % next character is a digit, add the character to the digit
                        % string and advance the position pointer
                    elseif p<L && isstrprop(s(p+1),'digit')
                        digstr=[digstr s(p+1)]; 
                        p=p+1;
                    else
                        error('Exponent in units string is not a proper number!')
                    end %if
                    % Convert the digit string to a number and multiply it by
                    % the increment in case the previous operator was a "/".
                    % Then set the increment back to zero
                    fundUnit.exp(ind)=inc*str2double(digstr);
                    inc=0;
            end %switch
            % Advance to the next character
            p=p+1;
        end %if
    end %while
    
    % If there is a token remaining, find the index of the fundamental unit and
    % increment the exponent appropriately
    if ~isempty(token)
        ind=strcmp(token,fundUnit.name);
        % Check to be sure that the token is a fundamental unit
        if ~any(ind)==1 && ~isempty(token)
            error(['The string "' token '" is not a fundamental unit!'])
        end %if
        fundUnit.exp(ind)=fundUnit.exp(ind)+inc;
    end %if
end % if
end % function

function s=fundUnits2str(F)
% Converts a fundamental units structure to a units string.
%   A fundamental units structure has the following form:
%        name: {'m'  'kg'  's'  'A'  'K'  'mol'  'cd'}
%         exp: [-2 2 0 0 -4 0 0]
% where the numbers represent the exponents of the various fundamental
% units. This function converts the structure to a single string. Units
% with positive exponents are placed in the string first using a multiply
% (*) symbol. Negative exponents are placed in the sting next with a divide
% (/) symbol.
% USEAGE:
%       s=fundUnits2str(F)

% Initialize the output string
s='';
% Find the locations of the exponents in "F" that are greater than zero
posExpLoc=find(F.exp>0);
% Find the locations of the exponents in "F" that are less than zero
negExpLoc=find(F.exp<0);
% If there are no positive exponents but it has a negative location, write
% it as a reciprocal
if isempty(posExpLoc) && ~isempty(negExpLoc)
    s=[s '1'];
end
% Write all exponents that are positive with the multiply operator
for I=posExpLoc
    % If the exponent is a positive "1", write it without an exponent
    % symbol. Otherwise write out the exponent.
    if F.exp(I)==1
        str=F.name{I};
        % If it is the last symbol in the numerator, write it without a
        % trailing multiply symbol. Otherwise write it with a trailing
        % multiply.
        if I==posExpLoc(end)
            s=[s str];  %#ok<*AGROW>
        else
            s=[s str '*']; 
        end % if
    else
        % If the exponent is larger than "1", write the fundamental unit
        % with an exponent.
        str=F.name{I};
        % If it is the last symbol in the numerator, write it without a
        % trailing multiply symbol. Otherwise write it with a trailing
        % multiply.
        if I==posExpLoc(end)
            s=[s str '^' num2str(F.exp(I))]; 
        else
            s=[s str '^' num2str(F.exp(I)) '*']; 
        end %if
    end %if
end % for

% If there are negative exponents write a division operator
if ~isempty(negExpLoc)
      s=[s '/'];
end %if

% Write all exponents that are negative with a division operator
for I=negExpLoc
    if F.exp(I)==-1
        str=F.name{I};
        % If it is the last symbol in the numerator, write it without a
        % trailing divide symbol. Otherwise write it with a trailing
        % divide.
        if I==negExpLoc(end)
            s=[s str]; 
        else
        s=[s str '/'];
        end % if
    else
        % If the exponent is larger than "1", write the fundamental unit
        % with an exponent.
        str=F.name{I};
        % If it is the last symbol in the denominator, write it without a
        % trailing divide symbol. Otherwise write it with a trailing
        % divide. Since it is being represented as a division, change the
        % sign of the exponent.
        if I==negExpLoc(end)
            s=[s str '^' num2str(-F.exp(I))];
        else
        s=[s str '^' num2str(-F.exp(I)) '/'];
        end %if
    end %if  
end % for
end % function

function [value,newStr]=removePrefix(str)
% This function identifies a prefix, removes it from the string, and
% returns the associated value and the new string with the prefix removed.

% Create the list of prefix names
prefixList={'exa','peta','tera','giga','mega','kilo','deka','deci',...
    'nano','pico','atto','yatto','zetta','hecto','centi','milli',...
    'micro','femto','zepto','yacto'};
% Create the list of values associated with the prefix
prefixValue=[1e18,1e15,1e12,1e9,1e6,1e3,10,0.1,1e-9,1e-12,1e-18,...
    1e24,1e21,1e2,1e-2,1e-3,1e-6,1e-15,1e-21,1e-24];
% Determine the length of the input string
L=length(str);
% Any prefix is at least three characters long. Check for a prefix only if
% there are three or more characters in the string
if L>2
    % Find the location(s) in the prefix list that match the first three
    % characters in the input string
    location=cellfun(@(x) ~isempty(x),strfind(prefixList,str(1:3)));
    % If any of the locations match, check whether the full prefix matches
    if any(location)
        prefix=prefixList{location};
        value=prefixValue(location);
        prefixLength=length(prefix);
        % If the string is longer than the prefix length and all of the
        % characters in the prefix match, return the new string with the
        % prefix removed and retain the value for the prefix
        if L>prefixLength && strcmpi(prefix,str(1:prefixLength))
            newStr=str(prefixLength+1:end);
        else
            % If the string does not match, set the new string to the
            % original string and the value to 1
            newStr=str;
            value=1;
        end %if
    else
        % If no location matches just return the new string as the origianl
        % string ad set the value to 1
        newStr=str;
        value=1;
    end %if
    % If no prefix is found, set the value to 1 and return the original
    % string
else
    newStr=str;
    value=1;
end % if
end % function
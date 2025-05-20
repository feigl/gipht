
% https://www.mathworks.com/help/matlab/matlab_oop/create-a-simple-class.html

% To use the class:
% 
% Save the class definition in a .m file with the same name as the class.
% 
% Create an object of the class.
% 
% Access the properties to assign data.
% 
% Call methods to perform operation on the data.
% 
% Create Object
% Create an object of the class using the class name:

a = IntPhaClass
% a = 
% 
%   IntPhaClass with properties:
% 
%     Value: []
% Initially, the property value is empty.
% 
% Access Properties
% Assign a value to the Value property using the object variable and a dot before the property name:

a.Value = pi

% To return a property value, use dot notation without the assignment:
a.Value
% For information on class properties, see Property Syntax.
 
% Call Methods
% Call the rad2int8 method on object a:

ipha.rad2int8
% ans =
% 
%     1.0500
% Pass the object as the first argument to a method that takes multiple arguments, as in this call to the multiplyBy method:

multiplyBy(a,3)
% ans =
% 
%     3.1416
% You can also call a method using dot notation:

a.multiplyBy(3)
% Passing the object as an explicit argument is not necessary when using dot notation. The notation uses the object to the left of the dot.
% 
% For information on class methods, see Method Syntax.
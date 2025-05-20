
% https://www.mathworks.com/help/matlab/matlab_oop/create-a-simple-class.html
classdef BasicClass
   properties
      Value {mustBeNumeric}
   end
   methods
      function r = roundOff(obj)
         r = round([obj.Value],2);
      end
      function r = multiplyBy(obj,n)
         r = [obj.Value]*n;
      end
   end
end

%classdef IntPhaClass
classdef IntPhaClass
    % define a class for wrapped phase as a signed 1-byte integer
    % -128 DN == -pi
    % +127 DN == +pi

    % https://www.mathworks.com/help/matlab/matlab_oop/create-a-simple-class.html

   properties
      Value {mustBeNumeric}
   end
   methods
      function i = rad2int8(obj)
         i = int8(256*round([obj.Value]/pi,2));
      end
      function c = arc(obj,a,b)
          % given phase values a and b each in radians on [-pi,pi]
          % return arc(a,b) on [0 pi]
          % see eq 2.3.13 page 19, Mardia and Jupp [2000]
          c = pi - abs(pi-abs(a-b));
      end
      function w = wrap(obj,a)
          % given phase values a in radians on [-pi,pi]
          % return arc(a,b) on [0 pi]
          % θ = 1 anglecomplexSinφ−φ ̃,
          % see eq 13 Feigl, K. L., and C. H. Thurber (2009), A method for
          % modelling radar interferograms without phase unwrapping:
          % application to the M 5 Fawnskin, California earthquake of 1992
          % December 4, Geophysical Journal International, 176, 491-504.
          % http://dx.doi.org/10.1111/j.1365-246X.2008.03881.x
          w = angle(complex(sin(a),cos(a)));
      end
   end
end

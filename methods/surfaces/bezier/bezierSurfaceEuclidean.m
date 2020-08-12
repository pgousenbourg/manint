function y = bezierSurfaceEuclidean(varargin)
% Computes a Bezier surace on the Euclidean space.
%
% function y = bezierSurfaceEuclidean(b,t1,t2)
%   returns y, the value of the Euclidean Bezier curve at t1 and t2. 
%   The curve is driven by the points b (cell).
%
% Original author: 
% 	Pierre-Yves Gousenbourger, May. 13, 2019.
% Contributors: 
%   ***
% Change log:
% 	May. 13, 2019 (PYG) - Original file.

% Manint - Copyright (C) <2014-2020> <Université catholique de Louvain (UCL), Belgique>
%	
% List of the contributors to the development of Manint: see AUTHORS file.
% Description and complete License: see LICENSE file.
%
% This program (Manint) is free software: 
% you can redistribute it and/or modify it under the terms of the 
% GNU General Public License as published by the Free Software Foundation, 
% either version 3 of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program (see COPYING file).  If not, 
% see <http://www.gnu.org/licenses/>.

  % create the parser
  ip = inputParser();
  addRequired(ip,'b');
  addRequired(ip,'t1');
  addRequired(ip,'t2');

  % parse the parser as a structure.
  parse(ip, varargin{:});
  vars = ip.Results;
  
  % variables
  b  = vars.b;
  t1 = vars.t1;
  t2 = vars.t2;
  
  assert(t1 >= 0 && t1 <= 1,'t1 must be between 0 and 1');
  assert(t2 >= 0 && t2 <= 1,'t2 must be between 0 and 1');
  
  deg = size(b,1)-1;
  
  % Curve reconstruction
  w = weights_bern(deg,t1,t2);
  w = w(:);
  b = b(:);
  
  % Weighted mean of the control points
  y = zeros(size(b{1}));
  for i = 1:length(b)
    y = y + w(i).*b{i};
  end
  
end


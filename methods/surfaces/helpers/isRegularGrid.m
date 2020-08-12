function bool = isRegularGrid(rootCoords)
% Checks if the grid given in ROOTCOORDS is regular.
%
% function isRegularGrid(ROOTCOORDS)
%    Returns a logical number saying if ROOTCOORDS specifies a regular
%    grid or not. ROOTCOORDS must be a [px2] matrix.
%
% See also: makeDomain
%
% Original author: 
% 	Pierre-Yves Gousenbourger, Jan. 07, 2020.
% Contributors: 
%
% Change log:
% 	Jan. 07, 2020 (PYG) - First version.

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

  % defense
  assert(size(rootCoords,2) == 2,'isRegularGrid:dimCheck','rootCoords must have two columns');
  
  % boolean return. Each occurence must appear the same number of time.
  bool = (isscalar(occurInArray(rootCoords(:,1))) && isscalar(occurInArray(rootCoords(:,2))));
end

% Extracts the number of occurences that appear in the same vector.
function A = occurInArray(X)
  X = hist(X(:));
  X = hist(X(X>0));
  A = X(X>0);
end

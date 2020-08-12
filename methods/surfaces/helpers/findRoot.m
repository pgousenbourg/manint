function idx = findRoot(X,Y,rootCoords)
% Computes the rootPoint index associated with the position (X,Y).
%
% idx = findRoot(X,Y,ROOTCOORDS)
%    finds the index of the rootCoord at position (X,Y). If none is
%    found, then an empty vector is returned.
%
% Inputs:  X,Y (scalars), positions in the domain.
%          ROOTCOORDS, a matrix with two columns with the (X,Y)-
%           positions of the grid coordinates.
% 
% Outputs: idx, the index of the rootCoord, in ROOTCOORDS.
%
% See also: findDomain, refineDomain, findNeighbors
%
% Original author: 
% 	Pierre-Yves Gousenbourger, Jan. 09, 2020.
% Contributors: 
%
% Change log:
% 	Jan. 09, 2020 (PYG) - First version.
%   Jan. 13, 2020 (PYG) - Defense and proper function.

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
  assert(size(rootCoords,2) == 2,'findRoot:dimCheck','rootCoords should have only 2 columns');
  assert(isscalar(X) && isscalar(Y),'findRoot:XY','X and Y must be scalar');
  
  % parameters
  tol = 1e-8;
  
  % find the indexes
  idx = find(abs(rootCoords(:,1) - X) < tol);
  if ~isempty(idx)
    idx = idx(find(abs(rootCoords(idx,2) - Y) < tol));
  end
end

function neighbors = findNeighbors(i,domains,rootCoords)
% Finds the four neighbors to domain i, at the same parenting level.
%
% neighbors = findNeighbors(I,DOMAINS,ROOTCOORDS)
%    finds the four neighbors [south,east,west,north] of domain i.
%    Returns 0 if no neighbor is found (for instance, if the domain i
%    is at the border of the meta-domain).
%
% Inputs:  I (scalar), the line of the domain that has to be refined.
%          DOMAINS, a [px5]-matrix of domains referring to the indexes
%           of ROOTCOORDS. The first line is the parenting.
%          ROOTCOORDS, a matrix with two columns with the (X,Y)-
%           positions of the grid coordinates.
% 
% Outputs: neighbors, is a 4-vector ordered like 
%            [south, east,west,north].
%          Returns 0 if no neighbor is found.
%
% See also: findDomain, refineDomain, findRoot, computeLimits
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
  assert(i >= 1 && i <= size(domains,1),'findNeighbors:dimCheck','Domain index out of bound.')
  assert(size(rootCoords,2) == 2,'findNeighbors:dimCheck','rootCoords must have two columns.')
  assert(size(domains,2) == 5,'findNeighbors:dimCheck','domains must have five columns.')
  
  % parameters
  X = unique(rootCoords(:,1));
  Y = unique(rootCoords(:,2));
  % minimal shift to go outside the limits of the domain, but remaining
  % in the limits of the neighbor.
  shiftX  = min(abs(X(1:end-1) - X(2:end)))./4;
  shiftY  = min(abs(Y(1:end-1) - Y(2:end)))./4;
  % limits and centers of the domain [south,east,west,north] 
  limits  = computeLimits(i,domains,rootCoords);
  center  = [limits(2)+limits(3),limits(1)+limits(4)]./2;
  % computation of the neighbor (take care of the limits of the domain)
  testPoints = zeros(4,2);
  testPoints(1,:) = [center(1),limits(1)-shiftY];
  testPoints(2,:) = [limits(2)-shiftY,center(2)];
  testPoints(3,:) = [limits(3)+shiftX,center(2)];
  testPoints(4,:) = [center(1),limits(4)+shiftY];
  
  neighbors = zeros(4,1);
  for k = 1:4
    try
      neighbors(k) = findDomain(testPoints(k,1),testPoints(k,2),domains,rootCoords);
    catch ME
      if strcmp(ME.identifier,'findDomain:outOfRange') 
        neighbors(k) = 0;
      else
        error('Critical error while finding the neighbor of a domain');
      end
    end
  end
end

function limits = computeLimits(i,domains,rootCoords)
% Computes the (X,Y) limits of the domain i.
%
% limits = computeLimits(I,DOMAINS,ROOTCOORDS)
%    computes the limits [ymin,xmin,xmax,ymax] of the domain i.
%
% Inputs:  I (scalar), the line of the domain that has to be refined.
%          DOMAINS, a [px5]-matrix of domains referring to the indexes
%           of ROOTCOORDS. The first line is the parenting.
%          ROOTCOORDS, a matrix with two columns with the (X,Y)-
%           positions of the grid coordinates.
% 
% Outputs: Limits are stored in a 4-vector ordered like limit
%            [south (ymin), east (xmin), west (xmax), north (ymax)].
%
% See also: findDomain, refineDomain, findRoot, findNeighbors
%
% Original author: 
% 	Pierre-Yves Gousenbourger, Jan. 09, 2020.
% Contributors: 
%
% Change log:
% 	Jan. 09, 2020 (PYG) - First version.
%   Jan. 13, 2020 (PYG) - Defense and proper function.

% note:
% /!\ WITH THIS HELPER, ASSOCIATED WITH FINDROOT, I COULD SIMPLIFY
%     FINDDOMAIN (I think ^^)

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
  assert(i >= 1 && i <= size(domains,1),'computeLimits:i','Domain index out of bound.')
  assert(size(rootCoords,2) == 2,'computeLimits:dimCheck','rootCoords must have two columns.')
  assert(size(domains,2) == 5,'computeLimits:dimCheck','domains must have five columns.')
  
  % function
  limits    = ones(4,1);
  corners   = rootCoords(domains(i,2:5),:);
  limits(1) = max(corners(1:2,2));
  limits(2) = max(corners([1,3],1));
  limits(3) = min(corners([2,4],1));
  limits(4) = min(corners(3:4,2));
end

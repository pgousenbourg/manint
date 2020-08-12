function domain = makeDomain(rootCoords)
% creates the matrix of the regular grid, where each line refers to
% the four indexes of the rootCoords vector of the four corners of the
% domain.
%
% function DOMAIN = makeDomain(ROOTCOORDS)
%    returns DOMAIN, the matrix of domains based on ROOTCOORDS.
%
% Inputs:  ROOTCOORDS must encapsulate the coordinates of a regular 
%          grid, given in a [(m*n)x2] matrix. Each line is a (X,Y)-
%          coordinate. Order has no importance; m and n correspond to
%          the number of X- and Y- coordinates in the regular grid.
%
% Outputs: DOMAIN is a [(m-1)*(n-1) x 5] matrix of the following form:
%          DOMAIN(:,1) is the parent column. It is set to zeros in this
%            initialization step. Parenting can is used when the grid is
%            refined (see also: refineGrid).
%          DOMAIN(:,2:5) contains the four indexes of the rootCoords
%            of the four corners of the domain.
%
% Corners of the domains are numbered in the following order 
% (grid-points are number as xy couples, hereunder). The same apprach is
% used for screening the domains.
%
%         01  ---  11       ->     w{3}  ---  w{4}
%          |       |                |          |
%          |       |                |          |
%         00  ---  10              w{1}  ---  w{2}
%
% See also: refineGrid, findDomain.
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

  % Defense on rootCoords
  assert(size(rootCoords,2) == 2,'surfaceFitting:dimCheck','rootCoords must have two columns');
  assert(isRegularGrid(rootCoords),'surfaceFitting:regularGrid','rootCoords must define a regular grid domain.');
  
  % Actual work
  [rC,I1] = sortCoords(rootCoords);
  [~,I2]  = sortCoords(rC,2); clear rC;
  
  % index swapping (useful for finding back the index in rootCoords)
  I = I1(I2);
  
  % number of coords in X and Y
  m = length(unique(rootCoords(:,1)));
  n = length(unique(rootCoords(:,2)));
  
  % fill in the domains
  domains = zeros((m-1)*(n-1),5);
  for j = 1:n-1
  for i = 1:m-1
    k   = (j-1)*(m-1) + i;  % index of the patch
    pos = (j-1)*m + i;      % first index taken in rC
    idx = [pos, pos+1, pos+m, pos+m+1]; % takes the indexes from rC
    domain(k,:) = [0 I(idx)']; % stores the indexes of dataCoords
  end
  end
  
end

function tangentValues = pointsLifting(M,rootPoints,points)
% Lifts the points in the tangent space of the last element of rootPoints 
% by transporting the points through all elements of rootPoints.
%
% function tangentValues = pointsLifting(M,rootPoints,points)
%    returns tangentValues, the representation in T_{rootPoints{end}}M 
%    of the points, lifted thanks to the M.log operator of M, a manifold
%    manopt-structure, and transported through M.transp.
%    The points are first lifted in the first element of rootPoints,
%    and then the tangent values are transported from element to element
%    of rootPoints, and the distance between each rootPoint is
%    accumulated at each transport.
%
% inputs:  M, a manopt manifold-structure containing the M.log operator.
%          rootPoints, a m-cell of point on the manifold M.
%          points, [DIM]-cell of points to be lifted to the tangent 
%            space at rootPoint.
%
% outputs: tangentValues, a [DIM]-cell containing the lifted version
%            the points to the tangent spaces at root.
%
% Original author: 
% 	Pierre-Yves Gousenbourger, Jan. 22, 2020.
% Contributors: 
%
% Change log:
% 	Jan. 22, 2020 (PYG) - First version.

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
  if ~iscell(rootPoints)
    temp = rootPoints;
    rootPoints = cell(1,1);
    rootPoints{1} = temp;
    clear temp;
  end
  if ~iscell(points)
    temp = points;
    points = cell(1,1);
    points{1} = temp;
    clear temp;
  end
  
  %rootPoints
  assert(isvector(rootPoints),'surfaceFitting:vectorCell','rootPoints must be a scalar or vectorial cell');
  
  % parameter
  sPoints     = size(points);
  sRootPoints = size(rootPoints);
  if sRootPoints > 1
    if isfield(M,'isotransp')
      tp = M.isotransp;
    else
      tp = M.transp;
    end
  end
  
  % method
  points = points(:);
  tangentValues = cell(size(points));
  for j = 1:length(points)
    temp = M.log(rootPoints{1},points{j});
    for k = 2:length(rootPoints)
      temp = tp(rootPoints{k-1},rootPoints{k},temp) + M.dist(rootPoints{k-1},rootPoints{k});
    end
    tangentValues{j} = temp;
    clear temp;
  end
  tangentValues = reshape(tangentValues,sPoints);
end

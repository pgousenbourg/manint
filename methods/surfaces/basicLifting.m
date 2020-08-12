function tangentValues = basicLifting(M,rootPoints,dataPoints)
% Lifts the dataPoints in the tangent space of each rootPoint.
%
% function tangentValues = basicLifting(M,rootPoints,dataPoints)
%    returns tangentValues, the representation in T_{rootPoints{i}}M of
%    the dataPoints, lifted thanks to the M.log operator of M, a manifold
%    manopt-structure.
%
% inputs:  M, a manopt manifold-structure containing the M.log operator.
%          rootPoints, a vector-cell of P elements with the rootPoint on 
%            the manifold M.
%          dataPoints, a vector-cell of Q points to be lifted to the 
%            tangent space at rootPoints.
%
% outputs: tangentValues, a [PxQ]-cell containing the lifted
%            version of the points to the tangent spaces at rootPoints.
%
% Original author: 
% 	Pierre-Yves Gousenbourger, Dec. 19, 2019.
% Contributors: 
%
% Change log:
% 	Dec. 19, 2019 (PYG) - First version.
%   Jan. 27, 2020 (PYG) - Renewed version based on pointsLifting

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

  if ~iscell(dataPoints)
    temp = dataPoints;
    dataPoints = cell(1,1);
    dataPoints{1} = temp;
    clear temp;
  end
  if ~iscell(rootPoints)
    temp = rootPoints;
    rootPoints = cell(1,1);
    rootPoints{1} = temp;
    clear temp;
  end
  
  assert(iscell(dataPoints) && isvector(dataPoints),'basicLifting:vectorCell','dataPoints must be in a vector-cell');
  assert(iscell(rootPoints) && isvector(rootPoints),'basicLifting:vectorCell','rootPoints must be in a vector-cell');
  
  % parameter
  sRoots = length(rootPoints);
  sData  = length(dataPoints);
  
  % method
  dataPoints = dataPoints(:);
  tangentValues = cell(sRoots,sData);
  for i = 1:sRoots
  for j = 1:sData
    tangentValues(i,j) = pointsLifting(M,rootPoints{i},dataPoints{j});
  end
  end
end

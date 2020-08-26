function tangentValues = coarseLifting(M,rootPoints,dataPoints,domains,rootCoords,dataCoords)
% Lifts the dataPoints in the tangent space of each rootPoint by
% by (i) summarizing them in a polynomial, (ii) following the shortest 
% path in the domain and (ii) reconstructing them afterwards.
%
% function tangentValues = coarseLifting(M,rootPoints,dataPoints,domains,rootCoords,dataCoords,coarseDegree)
%    returns tangentValues, the representation in T_{rootPoints{i}}M 
%    of the dataPoints, lifted thanks to the M.log operator of M, a 
%    manifold manopt-structure.
%    The dataPoints are first lifted in the closest rootPoint, and then
%    and then the tangent values are transported via M.transp() from 
%    rootPoint to rootPoint, following a shortest path in the domain.
%
% inputs:  M, a manopt manifold-structure containing the M.log operator.
%          rootPoints, a vector-cell of P elements with the rootPoint on 
%            the manifold M.
%          dataPoints, a vector-cell of Q points to be lifted to the 
%            tangent space at rootPoints.
%          domains, the domain where the rootPoints are set.
%          rootCoords, the coordinates of the rootPoints.
%          dataCoords, the coordinates of the dataPoints.
%
% outputs: tangentValues, a [PxQ]-cell containing the lifted
%            version of the points to the tangent spaces at rootPoints.
%
% Original author: 
% 	Pierre-Yves Gousenbourger, Jan. 20, 2020.
% Contributors: 
%
% Change log:
% 	Jan. 20, 2020 (PYG) - First version.
%   Jan. 27, 2020 (PYG) - Renewed version based on pointsLifting
%   Aug. 25, 2020 (PYG) - Solution based on a modified version of 
%                         transportLifting and pointsLifting

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

tangentValues = transportLifting(M,rootPoints,dataPoints,domains,rootCoords,dataCoords,1);

end

function tangentValues = transportLifting(M,rootPoints,dataPoints,domains,rootCoords,dataCoords,coarse)
% Lifts the dataPoints in the tangent space of each rootPoint by
% by following the shortest path in the domain.
%
% function tangentValues = transportLifting(M,rootPoints,dataPoints,domains,rootCoords,dataCoords)
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
%          coarse (optional) is a boolean that specifies if the points
%            must be approximated by a linear surface within a patch
%            (default: 0).
%            coarse and pointsCoords must be provided together.
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
%   Aug. 25, 2020 (PYG) - Adaptation to coarseLifting.

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

  if nargin < 7
    coarse = false;
  end
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
  
  % defense
  assert(iscell(dataPoints) && isvector(dataPoints),'basicLifting:vectorCell','dataPoints must be in a vector-cell');
  assert(iscell(rootPoints) && isvector(rootPoints),'basicLifting:vectorCell','rootPoints must be in a vector-cell');
  assert(size(domains,2) == 5,'transportLifting:dimCheck','Domains must have 5 columns.');
  assert(size(dataCoords,2) == 2,'transportLifting:dimCheck','dataCoords must have 2 columns.');
  assert(size(dataCoords,1) == length(dataPoints),'transportLifting:dimCheck','dataCoords and dataPoints must have the same length.');
  % defense for rootCoords (2 columns && size(1) == sRoots && max(domains) == sRoots)
  assert(size(rootCoords,2) == 2,'transportLifting:dimCheck','rootCoords must have 2 columns.');
  assert(size(rootCoords,1) == length(rootPoints),'transportLifting:dimCheck','rootCoords and rootPoints must have the same length.');
  assert(size(rootCoords,1) == max(max(domains(:,2:5))),'transportLifting:domains','domains refer to unexistant rootPoints.');
  
  
  % parameter
  sRoots = length(rootPoints);
  sData  = length(dataPoints);
  A      = domain2adjacency(domains,rootCoords);
  dataDomains = findDomain(dataCoords(:,1),dataCoords(:,2),domains,rootCoords);
  
  % stores the indexes of the routes in a cell
  routes = cell(sRoots,sRoots);
  tangentValues = cell(sRoots,sData);
  sDomains = unique(dataDomains); % domains on which there are dataPoints
  
  for out = 1:sRoots
    %for j = 1:sData
      %% step 1: detect the domain where the dataPoint lies
      %i = dataDomains(j);
      %% step 1.5: group dataPoints by domain (storage of indices);
      
    %end
    for j = 1:length(sDomains);
      % step 1: find the dataPoints who are in the domain
      i = sDomains(j);
      ptsIdx = find(dataDomains == i);
      points = dataPoints(ptsIdx);
      pointsCoords = dataCoords(ptsIdx,:);
      
      % step 2: detect the index of optimal first rootPoint
      corners = domains(i,2:5);
      cornerCoords = rootCoords(corners,:);
      diff    = sum((repmat(rootCoords(out,:),4,1) - cornerCoords).^2,2);
      [~,idx] = min(diff);
      in      = corners(idx);
      
      % step 3: compute the optimal route if necessary (i.e. if the route wasn't already computed)
      if isempty(routes{in,out})
        [~,routes{in,out}] = dijkstra(A,in,out);
      end
      nodes = rootPoints(routes{in,out});
      
      % step 4: transport and lifting
      tangentValues(out,ptsIdx) = pointsLifting(M,nodes,points,coarse,pointsCoords);
    end
  end
end

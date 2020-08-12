function [domains,rootCoords] = refineDomain(i,domains,rootCoords)
% Refines the domain i.
%
% [domains,rootCoords] = refineDomain(I,DOMAINS,ROOTCOORDS)
%    returns domains and rootCoords, where the line I of DOMAINS
%    is now parent of four subdomains (or only two, if DOMAINS(I,:) was
%    a semi-patch), and where the coordinates of the new grid
%    intersection are added to ROOTCOORDS.
%
% Inputs:  I (scalar), the line of the domain that has to be refined.
%           If I is a 2-vector, then it is considered as the position
%           (X,Y) within the domain that must be refined.
%          DOMAINS, a [px5]-matrix of domains referring to the indexes
%           of ROOTCOORDS. The first line is the parenting.
%          ROOTCOORDS, a matrix with two columns with the (X,Y)-
%           positions of the grid coordinates.
% 
% Outputs: domains, DOMAIN where DOMAIN(I,:) is now parent of subdomains
%            and where the neighbors are parents of semi-patches.
%          rootCoords, the ROOTCOORDS with new grid points.
%
% See also: findDomain, isSemi
%
% Original author: 
% 	Pierre-Yves Gousenbourger, Jan. 09, 2020.
% Contributors: 
%
% Change log:
% 	Jan. 09, 2020 (PYG) - First version.
%   Jan. 15, 2020 (PYG) - Merging makeNewDomains and completeSemiDomains
  
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
  assert(size(rootCoords,2) == 2,'refineDomain:dimCheck','dataCoords must have two columns');
  assert(size(domains,2) == 5,'refineDomain:dimCheck','domains must have five columns');
  if isscalar(i)
    assert(i > 0 && i < size(domains,1),'refineDomain:outOfRange','There is no line i in domains');
  end
  
  % Preprocessing: 
  % Transform i to a line-value if given as (x,y)-value.
  if length(i) == 2
    i = findDomain(i(1),i(2),domains,rootCoords);
  end
  
  % Parameters of the function
  tol = 1e-8;
  positions = {'SOUTH','WEST','EAST','NORTH'};
  
  if isSemi(i,domains,rootCoords)
    % Complete the semi-domains
    [domains,rootCoords] = makeNewDomains(i,domains,rootCoords);
    i = domains(i,1);
  elseif ~isParent(i,domains) && i~=0 % nothing to do for parent or level 0
    neighbors = findNeighbors(i,domains,rootCoords);
    for k = 1:length(neighbors)
      if neighbors(k) ~=0 % nothing to do if no neighbor
        if isSemi(neighbors(k),domains,rootCoords) 
          % refine the parent domain of semi-neighbors
          [domains,rootCoords] = refineDomain(neighbors(k),domains,rootCoords);
        end
      end
    end
    
    % Create four new domains for the patch i, if necessary
    [domains,rootCoords] = makeNewDomains(i,domains,rootCoords);
  end
  
  % Now, deal with the neighbors (creating semis if needed)
  neighbors = findNeighbors(i,domains,rootCoords);
  
  % For each neighbor:
  %   - if parent of four patches, then do nothing
  %   - if parent of two semis, then refine the domain
  %   - if normal, create semis
  
  % parameters of current patch
  limits    = computeLimits(i,domains,rootCoords);
  width     = abs(limits(2)-limits(3));
  height    = abs(limits(4)-limits(1));
  
  for k = 1:4
    if neighbors(k) ~= 0
      if isSemi(neighbors(k),domains,rootCoords)
        [domains,rootCoords] = refineDomain(neighbors(k),domains,rootCoords);
      elseif ~isParent(neighbors(k),domains)
        % parameters of the neighbor
        limNeigh = computeLimits(neighbors(k),domains,rootCoords);
        w = abs(limNeigh(2)-limNeigh(3));
        h = abs(limNeigh(4)-limNeigh(1));
        
        % do nothing if on smaller level
        switch positions{k}
        case {'SOUTH','NORTH'}
          sameLevel = abs(w - width) < tol;
        case {'WEST','EAST'}
          sameLevel = abs(h - height) < tol;
        end
        
        % apply semi patches if neighbor is on the same level
        if sameLevel
          [domains,rootCoords] = makeSemiDomains(neighbors(k),domains,rootCoords,positions{k});
        end
      end
    end
  end
end





% makeNewDomains
% Function used to complete a semi-domain or make four new domains
function [domains,rootCoords] = makeNewDomains(i,domains,rootCoords)
  if i~=0 && ~isParent(i,domains) % nothing to do for parents
    
    semi = isSemi(i,domains,rootCoords);
    if semi % detect the parent of semiPatches
      i = domains(i,1);
    end
    
    % parameters
    limits  = computeLimits(i,domains,rootCoords);
    center  = [limits(2)+limits(3),limits(1)+limits(4)]./2;
    corners = domains(i,2:5);
    
    % Creation of the possible new coordinates
    newCoords = [center(1),limits(1);
                 limits(2),center(2);
                 center(1),center(2);
                 limits(3),center(2);
                 center(1),limits(4)];
    
    % Storage of rootIdx among existing coords and new ones
    rootIdx = zeros(5,1);
    for k = 1:5
      idx = findRoot(newCoords(k,1),newCoords(k,2),rootCoords);
      % if no match, create a new coordinate.
      if isempty(idx)
        rootCoords = [rootCoords;newCoords(k,:)];
        idx = size(rootCoords,1);
      end
      rootIdx(k) = idx;
    end
    
    % Creation of new domains
    newDomains = [i corners(1) rootIdx(1) rootIdx(2) rootIdx(3);
                  i rootIdx(1) corners(2) rootIdx(3) rootIdx(4);
                  i rootIdx(2) rootIdx(3) corners(3) rootIdx(5);
                  i rootIdx(3) rootIdx(4) rootIdx(5) corners(4)];
    
    % Replace existing domains and insertion of new ones in the domains
    if semi
      domsIdx = find(domains(:,1) == i);
      domains(domsIdx,:) = newDomains(1:2,:);
      domains = [domains;newDomains(3:4,:)];
    else
      domains = [domains;newDomains];
    end
  end
end





% makeSemiDomains
% Function to create two brand new semiDomains in place of parent-domain i.
function [domains,rootCoords] = makeSemiDomains(i,domains,rootCoords,position)
  if i ~= 0 && ~isSemi(i,domains,rootCoords)
    limits  = computeLimits(i,domains,rootCoords);
    center  = [limits(2)+limits(3),limits(1)+limits(4)]./2;
    % the additional rootCoord depends on the position
    % the index-entry to be changed also depends on the position
    switch position
      case 'SOUTH'
        XY = [center(1),limits(4)];
        entry = [4 3];
      case 'WEST'
        XY = [limits(3),center(2)];
        entry = [4 2];
      case 'EAST'
        XY = [limits(2),center(2)];
        entry = [3 1];
      case 'NORTH'
        XY = [center(1),limits(1)];
        entry = [2 1];
      otherwise
        error('Unknown position for a patch');
    end
        
    % find the index of the additional rootCoord
    idx = findRoot(XY(1),XY(2),rootCoords);
    assert(~isempty(idx),'makeSemiDomains:emptyIndex','Fatal error: Impossible to find the additional rootCoord.');
    
    % Create new domains
    newDomains = repmat([i,domains(i,2:5)],2,1);
    % replace the correct entries
    newDomains(1,1+entry(1)) = idx;
    newDomains(2,1+entry(2)) = idx;
    % update domains
    domains = [domains;newDomains];
  end
end

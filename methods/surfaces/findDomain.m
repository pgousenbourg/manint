function idx = findDomain(X,Y,domains,rootCoords)
% finds the domain containing (X,Y).
%
% function IDX = findDomain(X,Y,DOMAINS,ROOTCOORDS)
%    returns IXD, the line of the domain from DOMAINS that contains the
%    coordinate (X,Y). The DOMAINS are a list of coordinates of 
%    ROOTCOORDS.
%
% Inputs: X,Y, two matrices of same size (DIM), with X and Y coordinates 
%           values.
%           If X and Y do not fit in the range prescribed by the
%           ROOTCOORDS, then 0 is returned instead, with a warning.
%         DOMAINS, the domain associated to the dataCoords.
%         ROOTCOORDS, the data coordinates of the corners of each patch
%           in the domain. There is no easy way to check that DOMAINS
%           and ROOTCOORDS are coherent with each other. So please: try
%           to use this function with care ;-).
%
% Output: IXD, a list of indexes of size DIM corresponding to the domains
%         containing (X,Y).
%
% See also: makeDomain, refineDomain.
%
% Original author: 
% 	Pierre-Yves Gousenbourger, Jan. 08, 2020.
% Contributors: 
%
% Change log:
% 	Jan. 08, 2020 (PYG) - First version.

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
  assert(size(rootCoords,2) == 2,'findDomain:dimCheck','dataCoords must have two columns');
  assert(size(domains,2) == 5,'findDomain:dimCheck','domains must have five columns');
  assert(length(size(X)) == length(size(Y)),'findDomain:dimCheck','X and Y must be the same DIM');
  
  % parameters
  sX = size(X);
  X = X(:); Y = Y(:);
  p = length(X);
  
  % defense (end)
  assert(length(X) == length(Y),'findDomain:dimCheck','X and Y must be of same size');
  assert(min(X) >= min(rootCoords(:,1)) && ...
         max(X) <= max(rootCoords(:,1)) && ...
         min(Y) >= min(rootCoords(:,2)) && ...
         max(Y) <= max(rootCoords(:,2)),'findDomain:outOfRange','coordinate (X,Y) out of the range of the domain');
  
  % allocation
  idx = zeros(p,1);
  
  % find domains for each pair (x,y)
  for i = 1:p
      % Search in the domains with smallest parenting level
      % output : variable "line"
      x = X(i);
      y = Y(i);
      line = findLocalDomain(x,y,domains,rootCoords,min(domains(:,1)));
      while isParent(line,domains)
        % do the search among the children.
        % output : update variable "line"
        line = findLocalDomain(x,y,domains,rootCoords,line);
      end
      
      % storage
      idx(i) = line;
  end
  
  idx = reshape(idx,sX);
end


% auxiliary function to search within domains of same level
%   - X,Y must be scalars
%   - domains is the full domain
%   - rootCoords are the full rootCoords
%   - parentLine is a scalar
function line = findLocalDomain(X,Y,domains,rootCoords,parentLine)
  % detect the local domains
  locDomainsIdx = find(domains(:,1) == parentLine);
  locDomains    = domains(locDomainsIdx,2:end); % no need to keep the parenting info
  
  if length(locDomainsIdx) == 2 % special guess for semi-patches
    % Compute the limits of one of the two semi-patches
    limits = computeLimits(locDomainsIdx(1),domains,rootCoords);
    % if not in the first, it is in the other
    if X >= limits(2) && X <= limits(3) && Y >= limits(1) && Y <= limits(4)
      line = locDomainsIdx(1);
    else
      line = locDomainsIdx(2);
    end
  else
    % detect the local roots within these local domains
    locRootsIdx   = unique(locDomains(:));
    locRoots      = rootCoords(locRootsIdx,:);
    
    % number of different rootPoints in X and Y
    m = length(unique(locRoots(:,1)));
    n = length(unique(locRoots(:,2)));
    
    % find the indexes with closest Y
    [rC,I] = sortCoords(locRoots,2);
    yIdx   = find(rC(:,2) - Y >= 0);
    yIdx   = yIdx(1:min(m,length(yIdx))); % FIND SOLUTION HERE !
                                          % m or m+1 ?
    
    % reduce the search space
    locRootsIdx = locRootsIdx(I(yIdx));
    locRoots    = rootCoords(locRootsIdx,:);
    
    % find the index with closest x
    [rC,I]  = sortCoords(locRoots,1);
    xIdx    = find(rC(:,1) - X >= 0);
    rootIdx = locRootsIdx(I(xIdx(1)));
    
    assert(isscalar(rootIdx),'rootIdx should be scalar');
    
    % Extract the line of the domain where rootIdx is the 4th corner point
    % (because it is the top right point)
    line = [];
    k = 4;
    while (isempty(line))
      line = locDomainsIdx(find(locDomains(:,k) == rootIdx));
      k = k-1;
    end
    
    % Decide when there is a draw. This happens when there is two semipatches
    %if length(line) == 2
      %[semi,pos] = isSemi(line,domains,rootCoords);
      %assert(strcmp(pos{1},pos{2}),'Fatal error: draw cannot happen for semi-patches at different positions. Contact the developper.');
      %rC = zeros(2,1);
      %switch pos{1}
        %case 'NORTH'
          %rC(1) = rootCoords(domains(line(1),1+2));
          %rC(2) = rootCoords(domains(line(2),1+2));
          %ref   = X;
        %case 'SOUTH'
          %rC(1) = rootCoords(domains(line(1),1+4));
          %rC(2) = rootCoords(domains(line(2),1+4));
          %ref   = X;
        %case 'WEST'
          %rC(1) = rootCoords(domains(line(1),1+4),2);
          %rC(2) = rootCoords(domains(line(2),1+4),2);
          %ref   = Y;
        %case 'EAST'
          %rC(1) = rootCoords(domains(line(1),1+3),2);
          %rC(2) = rootCoords(domains(line(2),1+3),2);
          %ref   = Y;
        %otherwise
          %error('Fatal error: unexpected output')
      %end
      %I = find(rC - ref >= 0);
      %line = line(I(1));
    %elseif length(line) > 2
      %error('Fatal error: there should never be a draw for more than two lines...');
    %end
  end
end

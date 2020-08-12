function [bool,varargout] = isSemi(i,domain,rootCoords)
% Checks if the ith line of domain is a semi-patch.
%
% BOOL = isSemi(I,DOMAIN,ROOTCOORDS)
%    Returns a logical number saying if the I^th line of DOMAIN is a
%    semi-patch.
%    If I is a vector, then the function will return a vector of
%    booleans checking every entry of I.
%
% [BOOL,POS] = isSemi(I,DOMAIN,ROOTCOORDS)
%    Also returns the position POS of the semi-patch with respect to
%    the refined patch that created him.
%    For instance, if the refined patch that created him is below,
%    POS will be 'NORTH'.
%    Possible values for POS is 'NORTH','WEST','SOUTH' and 'EAST'.
%
% See also: refineDomain, findDomain, makeDomain
%
% Original author: 
% 	Pierre-Yves Gousenbourger, Jan. 09, 2020.
% Contributors: 
%
% Change log:
% 	Jan. 09, 2020 (PYG) - First version.

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
  assert(size(rootCoords,2) == 2,'isSemi:dimCheck','dataCoords must have two columns');
  assert(size(domain,2) == 5,'isSemi:dimCheck','domains must have five columns');
  assert(isempty(find(i <= 0)),'isSemi:i','i must be greater than zero');
  
  % parameters
  si = size(i);
  i = i(:);
  
  % The actual work is done here
  for k = 1:length(i)
    bool(k) = length(find(domain(:,1) == domain(i(k),1))) == 2;
    if bool(k)
      assert(~isParent(i(k),domain),'isSemi:notParent','Fatal error: a semi-patch cannot be a parent');
    end
  end
  bool = reshape(bool,si);
  
  % Fill the pos if requested
  if nargout > 1;
    pos = cell(size(i));
    for k = 1:length(i)
      if bool(k)
        corners = domain(i(k),2:5);
        distSouth = abs(rootCoords(corners(1),1) - rootCoords(corners(2),1));
        distNorth = abs(rootCoords(corners(3),1) - rootCoords(corners(4),1));
        distWest  = abs(rootCoords(corners(1),2) - rootCoords(corners(3),2));
        distEast  = abs(rootCoords(corners(2),2) - rootCoords(corners(4),2));
        if distSouth < distNorth && distWest == distEast
          pos{k} = 'NORTH';
        elseif distSouth > distNorth && distWest == distEast
          pos{k} = 'SOUTH';
        elseif distSouth == distNorth && distWest < distEast
          pos{k} = 'EAST';
        elseif distSouth == distNorth && distWest > distEast
          pos{k} = 'WEST';
        else
          pos{k} = 'unexpectedBehavior';
        end
      else
        pos{k} = 'notSemi';
      end
    end
    if length(i) == 1;
      varargout{1} = pos{1};
    else
      varargout{1} = reshape(pos,si);
    end
  end
end

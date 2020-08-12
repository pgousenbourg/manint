function bool = isParent(i,domain)
% Checks if the ith line of domain is a parent of other domains.
%
% function isParent(DOMAIN,I)
%    Returns a logical number saying if the I^th line of DOMAIN has
%    associated childrens in DOMAIN.
%    If I is a vector, then the function will return a vector of
%    booleans checking every entry of I.
%
% See also: refineDomain, findDomain, makeDomain
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
  assert(size(domain,2) == 5,'isParent:dimCheck','domains must have five columns');
  assert(isempty(find(i <= 0)),'isParent:i','i must be greater than zero');
  
  % parameters
  si = size(i);
  i = i(:);
  
  bool = zeros(size(i));
  for k = 1:length(i)
    bool(k) = ~isempty(find(domain(:,1) == i(k)));
  end
  bool = reshape(bool,si);
end

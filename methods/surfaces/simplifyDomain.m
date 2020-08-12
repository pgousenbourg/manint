function simpleDomains = simplifyDomain(domains)
% Simplifies the domain by removing parents.
%
% function simpleDomains = simplifyDomain(domains)
%    returns simpleDomains, a simplified version of domains where all
%    the parents are removed.
%
% Inputs:  domains is a [px5] matrix.
% Outputs: simpleDomains, a [qx5] matrix where there is no parent
%          anymore (every parent line is now 0).
%
% Original author: 
% 	Pierre-Yves Gousenbourger, Jan. 16, 2020.
% Contributors: 
%
% Change log:
% 	Jan. 16, 2020 (PYG) - First version.

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
  assert(size(domains,2) == 5,'simplifyDomain:dimCheck','Domains must have 5 columns');
  
  simpleDomains = [];
  for i = 1:size(domains,1)
    if ~isParent(i,domains)
      simpleDomains = [simpleDomains; 0 domains(i,2:5)];
    end
  end
end

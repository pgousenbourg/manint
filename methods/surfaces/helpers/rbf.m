function y = rbf(r)
% Returns the value of the Radial Basis Function, used for Thin Plate
% Splines, at ray r.
%
% function y = rbf(r)
%    returns y, the value of the RBF at ray r, given by
%    rbf(t) = | 1/(16pi) r^2 log r^2  if r > 0
%             | 0                     else

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

  assert(length(r(:)) == 1,'r must be a scalar');
  if r > 0
    y = (1./(16*pi))*(r^2)*log(r^2);
  else
    y = 0;
  end
end

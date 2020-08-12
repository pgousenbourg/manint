function y = karcher(M,x,b,w)
% Computes the Riemannian center of mass of the points b weighted by w.
%
% function y = karcher(M,x,b,w)
%    returns y, the Karcher mean given by
%            y = sum_i w_i dist(b_i,x),
%    where dist is the geodesic distance. b is a cell of the same
%    dimension of w.
%
% Original author: 
% 	Pierre-Yves Gousenbourger, Nov. 07, 2018.
% Contributors: 
%   ***
% Change log:
% 	Nov. 07, 2018 (PYG) - Original file.

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

	b = b(:);
	w = w(:);
	
	elems = length(b);
	assert(elems == length(w),'The dimension of b and w must agree.');
	
	y = 0;
		
	for i = 1:elems;
    %size(x)
    %size(b{i})
		y = y + w(i)*M.dist(x,b{i}).^2;
	end
end

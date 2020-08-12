%% WEIGHTS_BERN_C1 -
%       Computation of weights based on the Bernstein representation
% 		for C^1 curves.

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

function b = weights_bern_c1(patch,deg,t,maxPatch)
	b = zeros(deg+1,1);
	
	
	% Extreme left (or down)
	if patch == 1
		b(1) = bern(deg,0,t);
		b(2) = bern(deg,1,t);
	else
		b(1) = 1/2*bern(deg,0,t);
		b(2) = bern(deg,1,t) + b(1);
	end
	
	% Inner points
	for i = 3:deg-2
		b(i) = bern(deg,i-1,t);
	end
	
	
	% Extreme right (or up)
	if patch == maxPatch
		b(deg)   = bern(deg,deg-1,t);
		b(deg+1) = bern(deg,deg,t);
	else
		b(deg+1) = 1/2*bern(deg,deg,t);
		b(deg)   = bern(deg,deg-1,t) + b(deg+1);
	end
end

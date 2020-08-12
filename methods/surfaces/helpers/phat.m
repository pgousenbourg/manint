% phat(p, pref, i, j): 
% 		 Compputes the k_th element of the P(p) vector.
% ------------------------------------------------------------
% This file is part of the project "C1 bezier paths on surfaces".
%
% INPUT: 	P 	 : [cell] the interpolation points of the problem.
% 					(1 x n)
% 
% 			k 	 : [int] The requested component of phat.
%
% OUTPUT: 	PHAT : [matrix] The k_th element of P(p)
% ------------------------------------------------------------
% Authors:
% 	Pierre-Yves Gousenbourger
% 	Paul Striewski
%
% Versions
% 	04/11/2015: First version.
% ------------------------------------------------------------

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

function phat = phat(p,k)
	
	% Useful data
	M = length(p); assert(M > 2, 'the length of p must be > 2');
	
	if M == 3
		phat = p{2} - p{1}/6 - p{3}/6;
	else
		if k == 2
			phat = 	p{2} - p{1}/6;
		elseif k == M-1
			phat =	p{k} - p{k+1}/6;
		else
			phat = 	p{k};
		end
	end
end

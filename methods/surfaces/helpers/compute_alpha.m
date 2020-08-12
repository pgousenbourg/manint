% compute_alpha(p, A, m, d): 
% 		 Computes the intermediate alpha of a C^2-Bezier spline in R^n,
% 		 such that 
% 				alpha = b_i^m = 'sum_{k=m-d}^{m+d} A_{k,m} p_k
% ------------------------------------------------------------
% This file is part of the project "C1 bezier paths on surfaces".
%
% INPUT: 	p : [cell] a line cell of points on a manifold M.
% 				(1 x M)
%
% 			A : [matrix] a matrix defining the linear combination.
% 				(M-2 x M-2)
%
%			m : [int] the position in the row
%
% 			d : [int] the size of the window of coefficients taken in A.
% 				(optional)
%
% OUTPUT: 	alpha : [cell] intermediate alphas in the computation of b.
%
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

function alpha = compute_alpha(p, A, m, d)
	
	% Parameters
	M 	  = length(p);
	
	% Computation
	switch m
	case 1
		alpha = p{1};
	case M
		alpha = p{M};
	otherwise
		upperBound = min(M-1,m+d);
		lowerBound = max(2,m-d);
		alpha = zeros(size(p{m}));
		for k = lowerBound:upperBound
			alpha = alpha + A(k-1,m-1)*phat(p,k); 
		end
	end
end

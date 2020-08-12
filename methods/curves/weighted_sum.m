function x = weighted_sum(w,y,M,u)
% Computes a weighted sum.
% 
% function x = weighted_sum(w,y)
% Computes the weighted sum x = sum_i ( w_i*y_i ). More specifically, it
% solves the problem
%      x = argmin_z sum_i ( w_i * dist(z,y_i) )
%
% Input: w, a vector of weights;
%        y, the data points in a m-by-n-by-p, where [m,n] is the
%           the dimension of the space and p is the number of points.
%           length(w) must be the same size as p.
%        M, the manifold structure containing differential geometry.
% 		 u, the reference point as m-by-n matrix. If not given, it will
% 			be the first point.
%
% Output: x, the summary vector of size [m,n].
%
% Original author: 
%   Estelle Massart, Oct. **, 2016
% Contributors: 
%	Pierre-Yves Gousenbourger, Oct. 17, 2016
%
% Change log:
% 	Oct. 17, 2016 (PYG):
%      Preparation of the document for integration into the framework
%      for approximation methods with Bezier functions.
%   Oct. **, 2016 (EM):
%      First version.

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
	
	if nargin < 3
		error('Not enough input arguments');
	else if nargin < 4
		u = y(:,:,1);		
	end
	
	% TO_DO_ESTELLE_MASSART
	x = zeros(size(y(:,:,1)));
	
	% Project on the tangent space of u (if defined)
	for i = 1:size(y,3)
		y(:,:,i) = M.log(u,y(:,:,i));
	end
	
	% Weighted mean
	for i = 1:length(w)
		x = x + w(i)*y(:,:,i);
	end
	
	% Back projection
	x = M.exp(u,x);
	
end

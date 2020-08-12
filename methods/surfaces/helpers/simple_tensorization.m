function b = simple_tensorization(M, dataPoints, A, threshold)
% Computes the control points of a C^2-Bezier spline in R^n.
%
% function b = simple_tensorization(M,dataPoints,A,threshold)
%    returns b, a cell containing data points of a C2 Bezier spline in
%    the Euclidean space. The manifold M is used to verify the distance
%    of the points wrt the geodesic. A is the matrix that generates the
%    control points. Threshold limits the input from the matric.
%
% Original author: 
% 	Pierre-Yves Gousenbourger, Nov. 02, 2015.
% Contributors: 
%	Paul Striewski.
% Change log:
% 	Nov. 05, 2018 (PYG) - Integration to the framework.
% 	Nov. 07, 2018 (PYG) - Integration of the Manopt framework.

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

	m = length(dataPoints);

	b = cell(3*m-2,1);

	for i = 1:m-1
	for k = 0:3
		b{3*(i-1)+k+1} = bki(M,dataPoints,A,k,i,threshold);
	end
	end

	% Sanity checks - Interpolation conditions
	for i=4:3:3*m-3
		av = M.pairmean(b{i-1}, b{i+1}); 
		dist = acos(dot(av,b{i}));
		if dist > 1e3
			warning('The distance between average and cp is high: %d\n', dist);
		end
	end
end

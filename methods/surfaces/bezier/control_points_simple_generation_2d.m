function pb = control_points_simple_generation_2d(pb)
% Returns the control points of an interpolating Bezier surface to
% a set of data points pb.dataPoints(i,j) at integer times 
% (t1,t2) := (i-1,j-1). On the Euclidean space, the obtained surface is
% C^2.
%
% function pb = control_points_simple_generation_2d(pb)
%    returns pb, the structure of the problem, updated with the control
% 	 points pb.controlPoints (in a 2d-cell).
%
%    The control points are obtained in a way that the mean squared 
%    acceleration of the surface is minimized on the Euclidean
%    space, and such that the surface can be computed in a C^2 manner.
%
% Input: pb, the structure of the interpolating problem, containing
% 			pb.dataPoints (2-cell with data points)
% 			pb.degree 	  (the degree of the curve: default 3)
%           pb.threshold  (the threshold on the weights of the final
% 	                       combination)
%           pb.patchM     (the number of patches in the first coordinate
%                           - rows, counted from top to bottom)
% 			pb.patchN     (the number of patches in the second coordinate
% 						    - columns, counted from left to right)
% 			pb.manifold	  (the manifold structure from Manopt)
%           
% Output: pb, the structure of the approximation problem, edited
%         with the field pb.controlPoints, the (deg+1)^2 control points 
%         of each patch, stored in a (M*deg+1, N*deg+1)-cell).
%
% Original author: 
% 	Pierre-Yves Gousenbourger, Nov. 04, 2015.
% Contributors: 
%	Paul Striewski.
% Change log:
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
	
	% defense
	if ~isfield(pb,'degree'); pb.degree = 3; end
	if ~isfield(pb,'patchM'); pb.patchM = size(pb.dataPoints,1)-1; end;
	if ~isfield(pb,'patchN'); pb.patchN = size(pb.dataPoints,2)-1; end;
   
	% parameters
  dataPoints  = pb.dataPoints;
  M           = pb.M;
  d           = pb.degree;
  [m,n]   	  = size(dataPoints);
  [r,c] 	    = size(dataPoints{1,1});
	if ~isfield(pb,'threshold'); pb.threshold = max(m,n); end;
	threshold   = pb.threshold; % the number of elements around the diag of Am, An.
	
	% preallocation
	b = cell(pb.patchM*d + 1,pb.patchN*d + 1);
	
	% defense
	assert(m > 2, 'There must be at least 3 interpolation points in the x-direction.');
	assert(n > 2, 'There must be at least 3 interpolation points in the y-direction.');
	assert(d == 3, 'The degree of the curve must be 3 - no working code available for another d');
	
	
	
	% Helpers
	Am = inv((1/6)*(diag(4*ones(1,m-2),0) + diag(ones(1,m-3),1) + diag(ones(1,m-3),-1)));
	An = inv((1/6)*(diag(4*ones(1,n-2),0) + diag(ones(1,n-3),1) + diag(ones(1,n-3),-1)));

	% Preparation of structures
	%b 	 	= cell(d*m - 2, d*n - 2);
	
	
	
	% ==== Computation in 2d ===========================================
	% Interior points
	for i = 1:m-1
	for j = 1:n-1
		for k = 0:d
		for l = 0:d
			b{3*(i-1) + k + 1, 3*(j-1) + l + 1} = bklij(M,dataPoints, Am, An, k, l, i, j, threshold);
		end
		end
	end
	end
	% ==== End computation in 2d =======================================
	
	
	%Interpolation points
	b(1:3:end,1:3:end) = dataPoints;
	
	% Control points at the interfaces
	% x direction
	for i = 1:m-2
		for l = 1:pb.patchN
			yshift = 3*l;
			b{3*i+1, yshift-1} = M.pairmean(b{3*i, yshift-1}, b{3*i + 2, yshift-1});
			b{3*i+1, yshift}   = M.pairmean(b{3*i, yshift  }, b{3*i + 2, yshift});
		end
	end
		
	% y direction
	for j = 1:n-2
		for k = 1:pb.patchM
			xshift = 3*k;
			b{xshift-1,3*j+1} = M.pairmean(b{xshift-1,3*j}, b{xshift-1,3*j+2});
			b{xshift  ,3*j+1} = M.pairmean(b{xshift  ,3*j}, b{xshift  ,3*j+2});
		end
	end	
		
	pb.controlPoints = b;

end

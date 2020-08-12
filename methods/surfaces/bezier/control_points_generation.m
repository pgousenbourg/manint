function pb = control_points_generation(pb)
% Returns the control points of an interpolating Bezier surface to
% a set of data points pb.dataPoints(i,j) at integer times 
% (t1,t2) := (i-1,j-1).
%
% function pb = control_points_generation(pb)
%    returns pb, the structure of the problem, updated with the control
% 	 points pb.controlPoints (in a 2d-cell).
%
%    The control points are obtained by a such that the following mean
%    squared acceleration of the surface is minimized on the Euclidean
%    space, and such that the surface can be computed in a C^1 manner.
%
% Input: pb, the structure of the interpolating problem, containing
% 			pb.dataPoints (2-cell with data points)
% 			pb.degree 	  (the degree of the curve: default 3)
% -- RM     pb.threshold  (the threshold on the weights of the final
% 	                       combination)
%       pb.patchM     (the number of patches in the first coordinate
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
% 	Benedikt Wirth, Jul. 14, 2015.
% Contributors: 
%	Pierre-Yves Gousenbourger, Jul. 14, 2015.
% Change log:
% 	Jul. 27, 2015 (PYG) - Running version for matrix manifolds.
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
  [m,n]   	  = size(pb.dataPoints);
  [r,c] 	    = size(dataPoints{1,1});
	
	% preallocation
	pb.controlPoints = cell(pb.patchM*d + 1,pb.patchN*d + 1);
	
	% parallel transport
  if ~isfield(M,'isotransp'); M.isotransp = M.transp; end
	parTp = @(i,j,k,l,v) M.isotransp(dataPoints{i,j},dataPoints{k,l},v);


	% Step 1 -- compute Z:
	% Logarithms of interpolation points to their neighbours
	% ----------------------------------------------------------
	Z = zeros(m,n,d,d,r,c);
	for i = 1:m
	for j = 1:n
		for k = 1:d
		for l = 1:d
			% Condition for the corners
			% if ( abs(2*(i+k-3)-n+1) <= n-1 ) && ( abs(2*(j+l-3)-m+1) <= m-1 )
			if (i*k > 1) && (i*k < 3*m) && (j*l > 1) && (j*l < 3*n)
				speed = M.log(dataPoints{i,j},dataPoints{i+k-2,j+l-2});
				Z(i,j,k,l,:,:) = speed;
			end
		end
		end
	end
	end


	% Step 2 -- Compute the alpha
	% The coefficients of the optimization problem
	% ----------------------------------------------------------
	if exist('BiStiffCoeff.mat','file') == 2
		% Load the coefs.
		fprintf('   Stiffness coefficients already exist... ')
		coeff = load('BiStiffCoeff');
		coeff = coeff.coeff;
	else
		% Coeffs for optimization and storage
		fprintf('   Stiffness coefficients must be computed... ')
		coeff = assembleBiStiffMat;
		save('BiStiffCoeff.mat','coeff');
	end
	disp('local stiffness matrix assembled.');



	% Step 3 -- Elastop.
	% From the Z values and the coefficients, computation of the
	% left and right members of the linear system in U.
	
	% "Minus" because Z is computed reversely
	b = 	- reshape( ...
				compute_b( Z, coeff, parTp ), ...
				[m*n*3*r*c 1] );		
				
	A = @(U)  reshape( ...
				compute_A( reshape( U, [m n d r c] ), coeff, parTp ), ...
				[m*n*d*r*c 1] );


	% Optimal velocity minimizing the energy on the Euclidean space
	U = pcg( A, b, 1e-6, 1000 );




	% compute control points as Exp y +/- v
	y = zeros(m,n,r,c);
	for i = 1:m
	for j = 1:n
		y(i,j,:,:) = dataPoints{i,j};
	end
	end

	U = reshape( U, [m n d r c] );
	v1 = reshape(U(:,:,1,:,:),[m n r c]);
	v2 = reshape(U(:,:,2,:,:),[m n r c]);
	v3 = reshape(U(:,:,3,:,:),[m n r c]);
	cp(:,:,2,2,:,:) = Expo(M,y,zeros(size(v1)));
	cp(:,:,3,2,:,:) = Expo(M,y,v1);
	cp(:,:,1,2,:,:) = Expo(M,y,- v1);
	cp(:,:,2,3,:,:) = Expo(M,y,v2);
	cp(:,:,2,1,:,:) = Expo(M,y,- v2);
	cp(:,:,3,3,:,:) = Expo(M,y,v3);
	cp(:,:,1,3,:,:) = Expo(M,y,- v3 + 2*v2);
	cp(:,:,3,1,:,:) = Expo(M,y,- v3 + 2*v1);
	cp(:,:,1,1,:,:) = Expo(M,y,v3 - 2*v2 - 2*v1);

	% Cp is made of [m,n,d,d,r,c] corresponding to:
	% 	m: index x of the dataPoints
	%   n: index y of the dataPoints
	% 	d,d: around the dataPoints
	% 	r,c: dimension of the dataPoints.

	% Reorder this into control points. There is too much control points
	% because fictive control points around the borders are also computed.
	bb = cell(m*d,n*d);
	for i=1:m;
	for j=1:n;
		
		bb{(i-1)*d+1,(j-1)*d+1} = reshape( cp(i,j,1,1,:,:),[r,c]);
		bb{(i-1)*d+2,(j-1)*d+1} = reshape( cp(i,j,2,1,:,:),[r,c]);
		bb{(i-1)*d+3,(j-1)*d+1} = reshape( cp(i,j,3,1,:,:),[r,c]);
		bb{(i-1)*d+1,(j-1)*d+2} = reshape( cp(i,j,1,2,:,:),[r,c]);
		bb{(i-1)*d+2,(j-1)*d+2} = reshape( cp(i,j,2,2,:,:),[r,c]);
		bb{(i-1)*d+3,(j-1)*d+2} = reshape( cp(i,j,3,2,:,:),[r,c]);
		bb{(i-1)*d+1,(j-1)*d+3} = reshape( cp(i,j,1,3,:,:),[r,c]);
		bb{(i-1)*d+2,(j-1)*d+3} = reshape( cp(i,j,2,3,:,:),[r,c]);
		bb{(i-1)*d+3,(j-1)*d+3} = reshape( cp(i,j,3,3,:,:),[r,c]);
		
	end
	end

	%% Structure update
	pb.controlPoints = bb(2:end-1,2:end-1);


end



function y = Expo(M,p,v)
	
	[m,n,r,c] = size(p);
	y = zeros(size(p));
	
	for i = 1:m
	for j = 1:n
		pp = reshape(p(i,j,:,:),[r,c]);
		vv = reshape(v(i,j,:,:),[r,c]);
		y(i,j,:,:) = M.exp(pp,vv);
	end
	end
end

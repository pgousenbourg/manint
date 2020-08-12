% FUNCTION COMPUTE_B(Z, ALPHA, PARTP): 
% 		 Computes the B member of the linear system to solve
% 		 when generating the control points.
% ------------------------------------------------------------
% This file is part of the project "C1 bezier paths on surfaces"
% and is intended in computing the velocity on the Bezier surface.
%
% INPUT: 	Z 	 : the distances between an interpolation point and
% 				   its 8 neighbours. Matrix(m,n,3,3,dim1,dim2).
%
% 				   Z(m,n,i,j,:,:) is interpreted as the distance between
% 				   the interpolation point p(m,n) and the neighbours 
% 				   p(m+j-2,n+i-2).
% 
% 			ALPHA: the Bistiff coefficients of the problem.
%
% 			PARTP: the parallel transport on the manifold.
%
% OUTPUT: 	RES: The result b stored as a matrix(m,n,3,dim1,dim2), 
% 				 representing the solution in b with respect to the 3 
% 				 control points to optimize.
%
% 				 The result res(m,n,i,:,:) is the result with respect to 
% 				 u(i), where u(i) is a control point to optimize, mapped
% 				 on the interpolation point p(m,n) respectively on:
% 				 i=1 - the east of p(m,n);
% 				 i=2 - the north of p(m,n);
% 				 i=3 - the north east of p(m,n).
% ------------------------------------------------------------
% Versions
% 	30/03/2015: first version.
% 	14/07/2015: update to matrix manifolds
% 	27/07/2015: Running version for matrix manifolds.
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

function res = compute_b( Z, alpha, parTp )

	[n,m,~,~,dim1,dim2] = size(Z);

	% initialisation
	energy 	= 0;
	auxRes 	= zeros(n,m,3,3,dim1,dim2);  	% energy derivative wrt Unear
	res 	= zeros(n,m,  3,dim1,dim2);     % energy derivative wrt U

	% iterate over all patches
	for i = 0:n-2		% Patch x
	for j = 0:m-2		% Patch y

		% iterate over the four corners of that patch
		for x = 1:2			% Corner x
		for y = 1:2			% Corner y

			
			% Operator Z:
			% -----------
			% This step is: move from U to V
			% (transport vectors to tangent space at that corner)
			% Here, we store in vLoc the shift Z needed in tildeT 
			% when U is moved from his interpolation point to the 
			% corner x,y of the patch n,m. This is done for all 16
			% control points of the patch, indeed.
			
			vLoc = zeros(4,4,1,1,dim1,dim2);
			for k = 1:4
			for l = 1:4
			
				pInd = round(([k,l]-1)/3)+1;  	% closest interpolation point
				UInd = [k,l]-1-3*(pInd-1)+2; 	% local vector index at that interpolation point (1,2 or 3)
				
				% Z(n,m,i,j) = from p(n,m) to p(n+i-2,n+j-2).
				vLoc(k,l,1,1,:,:) = Z(  i+x 		  ,	j+y, ...
										pInd(1) - x+2 ,	pInd(2) - y + 2,...
										:,:);
			end
			end
			
			
			
			% Operator L:
			% -----------
			% compute the biharmonic energy of the patch and linear 
			% operator, all based on tangent space at current corner.
			
			locRes = zeros(4,4,1,1,dim1,dim2);
			for k = 1:4
			for l = 1:4
			for o = 1:4
			for p = 1:4
				%energy = energy + biStiffCoeff(k,l,o,p) * dotProd(squeeze(vLoc(k,l,1,1,:)),squeeze(vLoc(o,p,1,1,:)));
				locRes(k,l,1,1,:,:) = locRes(k,l,1,1,:,:) + alpha(k,l,o,p) * vLoc(o,p,1,1,:,:);
				locRes(o,p,1,1,:,:) = locRes(o,p,1,1,:,:) + alpha(k,l,o,p) * vLoc(k,l,1,1,:,:);
			end
			end
			end
			end
			
			
			% Reversed operation: transport the vectors to their
			% original corner
			for k = 1:4
			for l = 1:4
				pInd = round(([k,l]-1)/3)+1;  % closest interpolation point
				UInd = [k,l]-1-3*(pInd-1)+2;  % local vector index at that interpolation point
				
				
				%sizeAux = size(auxRes(i+pInd(1),j+pInd(2),UInd(1),UInd(2),:,:));
				vector = reshape(locRes(k,l,1,1,:,:),[dim1,dim2]);
				locResTP = parTp(i+x,j+y,i+pInd(1),j+pInd(2),vector); 
				
				%locResTP = reshape( parTp(i+x,j+y,i+pInd(1),j+pInd(2),locRes(k,l,1,1,:,:)),sizeAux); 
				
				auxResTemp = reshape(auxRes(i+pInd(1),j+pInd(2),UInd(1),UInd(2),:,:),[dim1,dim2]) + locResTP;
				auxRes(i+pInd(1),j+pInd(2),UInd(1),UInd(2),:,:) = auxResTemp;
				%auxRes(i+pInd(1),j+pInd(2),UInd(1),UInd(2),:,:) = ...
						%auxRes(i+pInd(1),j+pInd(2),UInd(1),UInd(2),:,:) ...
					%+ 	locResTP;
				
			end
			end

		end		% Corner y
		end		% Corner x

	end		% Patch y
	end		% Patch x



	% compute derivative wrt U by chain rule
	for i = 1:n
	for j = 1:m
		res(i,j,1,:,:) = auxRes(i,j,3,2,:,:) - auxRes(i,j,1,2,:,:) + 2*auxRes(i,j,3,1,:,:) - 2*auxRes(i,j,1,1,:,:);
		res(i,j,2,:,:) = auxRes(i,j,2,3,:,:) - auxRes(i,j,2,1,:,:) + 2*auxRes(i,j,1,3,:,:) - 2*auxRes(i,j,1,1,:,:);
		res(i,j,3,:,:) = auxRes(i,j,3,3,:,:) - auxRes(i,j,1,3,:,:) -   auxRes(i,j,3,1,:,:) +   auxRes(i,j,1,1,:,:);
	end
	end

end

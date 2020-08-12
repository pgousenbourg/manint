function controlPoints = bezierSurfaceSplineGeneration(dataPoints,dataCoords,lambda,varargin)
% Returns the coefficients conducting the Euclidean Bezier Surface.
%
% function controlPoints = bezierSurfaceSplineGeneration(dataPoints,dataCoords)
%    returns controlPoints, a [MxN]-cell containing the control points
%    of the bidimentional Bezier Surface interpolating the dataPoints{i} 
%    at dataCoords(i). The Bezier surface is of degree 3.
%
% function controlPoints = bezierSurfaceSplineGeneration(dataPoints,dataCoords,lambda)
%    returns controlPoints, a [MxN]-cell containing the control points
%    of the bidimentional Bezier Surface interpolating the dataPoints{i} 
%    at dataCoords(i) with a fitting parameter lambda (the smaller, the closer).
%
% function controlPoints = bezierSurfaceSplineGeneration(dataPoints,dataCoords,lambda,'rootCoords',ROOTCOORDS)
%    Specifies the position of the coordinates of the corners of the patches.
%    By default, they are considered as integers from [0,patchM]x[0,patchN],
%    at each integer value.
%
% function controlPoints = bezierSurfaceSplineGeneration(dataPoints,dataCoords,lambda,'options',OPTIONS)
%    adds options to the generation. OPTIONS must be a structure. It can
%    contain the following fields
%       OPTIONS.patchM, the number of patches along the first direction;
%       OPTIONS.patchN, the number of patches along the second direction;
%       OPTIONS.thershold, the threshold for weighting means (default: 1e-6);

% inputs:  dataPoints is a n-dimensional cell. Stores the data points
%            of size [p,q].
%          dataCoords is a [nx2]-matrix. Stores the coordinates at which
%            the data points must be fitted.
%          lambda     is a scalar (optional). The smaller, the closer to
%            interpolation the curve will be.
%
% outputs: controlPoints is a [MxN]-cell with the values of the control
%          points.
%
% Original author: 
%   Benedikt Wirth, Jan. 31, 2018
% Contributors: 
%   Pierre-Yves Gousenbourger, Apr. 26, 2018
% Change log:
% 	Apr. 26, 2018 (PYG) - Correction of indexes in the construction of C
%   May. 06, 2020 (PYG) - reduction to Euclidean case and degree 3

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
  
  % Defense
  assert(iscell(dataPoints),'surfaceFitting:dataNotInCell','The data points must be stored in a cell');
  dataPoints = dataPoints(:);
  
  assert(size(dataCoords,2) == 2,'surfaceFitting:coordsNotIn2D','The dataCoords must be stored in a nx2 matrix');
  assert(size(dataCoords,1)==length(dataPoints),'surfaceFitting:dimCheck','There must be as many dataPoints as dataCoords');
  if nargin < 3
    lambda = 1e8;
  end
  
  if nargin >= 3
    % create the parser
    ip = inputParser();
    addOptional(ip,'options',[]); % empty by default
    addOptional(ip,'rootCoords',[]); % empty by default
    % parse the parser as a structure.
    parse(ip, varargin{:});
    vars = ip.Results;
  end
  
	% Data processing
  if (~isfield(vars.options,'patchM') && isfield(vars.options,'patchN')) || (isfield(vars.options,'patchM') && ~isfield(vars.options,'patchN'))
    error('Please provide patchM AND patchN in the option field, or none of them');
  end
  if ~isfield(vars.options,'patchM') && ~isfield(vars.options,'patchN');
    if ~isempty(vars.rootCoords) && isRegularGrid(vars.rootCoords)
      vars.options.patchM = length(unique(rootCoords(:,1))) - 1;
      vars.options.patchN = length(unique(rootCoords(:,2))) - 1;
    else
      warning('rootCoords are not a regular grid and no patch number is given. Default: 5. Ignoring the coordinates.');
      vars.options.patchM = 5;
      vars.options.patchN = 5;
      vars.rootCoords = [linspace(min(dataCoords(:,1)),max(dataCoords(:,1)),6)' linspace(min(dataCoords(:,2)),max(dataCoords(:,2)),6)'];
    end
  end
  rootCoords = vars.rootCoords;
  % reparameterizing the data points on [0,patchM]x[0,patchN]
  dataCoords = (dataCoords - min(rootCoords))./(max(rootCoords) - min(rootCoords)).*([vars.options.patchM vars.options.patchN]);
  if ~isfield(vars.options,'threshold')
    vars.options.threshold = 1e-6;
  end

  options = vars.options;


  % Parameters of the method
	M = vars.options.patchM;
	N = vars.options.patchN;
	deg = 3;
	thresh = vars.options.threshold;
	ndp = size(dataCoords,1); % number of data points
	ncp = M*N*(deg+1)^2; 	    % number of control points
	
	
	% ==================================================================  
	% Matrix A (stiff)
	% ----------------
	% Assemble the stiffness matrix A, for each patch, in terms of the
	% Bezier control points, from 1/2 * x^T A x.
	%
	% A_ijop = \int_{[0,1]^2} [B_j''B_i B_j'B_i';B_j'B_i' B_jB_i'']:[B_o''B_p B_o'B_p';B_o'B_p' B_oB_p''] dt1 dt2
	% ==================================================================
	aux1 = zeros(deg+1,deg+1);
	aux2 = aux1;
	aux3 = aux1;
	for i = 0:deg
	for j = i:deg
		% This was obtained with symbolic calculation. Bin is the Bernstein polynomial: 
		% In maple: 
		% 	Bin := nchoosek(n,i)*x^i*(1-x)^(n-i)
		
		% Part with no derivative:
		% int_[0,1] B_in (t) B_jn(t)
		%
		% In maple: 
		% 	int(Bin*subs(i=j,Bin),x=0..1);
		%aux1(i+1,j+1) = nchoosek(deg,j) * nchoosek(deg,i) / nchoosek(2*deg+1,i+j) / (2*deg+1-i-j); 
		aux1(i+1,j+1) = f_aux1(i,j,deg);
		
		% Part with first derivative
		% int_[0,1] B'_in(t) B'_jn(t)
		%
		% In maple: 
		% 	simplify(
		%    	int(diff(Bin,x)*subs(i=j,diff(Bin,x)),x=0..1) 
		%       assuming i::integer, n::integer, j::integer
		%   );
		%aux2(i+1,j+1) = nchoosek(deg,j)*nchoosek(deg,i)/nchoosek(2*deg-4,i+j-2)*(deg-1)*(deg*(i+j)-2*i*j-deg*(i-j)^2)/((2*deg-3)*(2*deg-2)*(2*deg-1)); 
		aux2(i+1,j+1) = f_aux2(i,j,deg);
		
		% Part with the second derivative.
		% int_[0,1] B''_in(t) B''_jn(t)
		% 
		% In maple: 
		% 	simplify(
		% 		int(diff(Bin,x$2)*subs(i=j,diff(Bin,x$2)),x=0..1) 
		% 		assuming i::integer, n::integer, j::integer
		% 	);
		%aux3(i+1,j+1) = ((n-3)*((i^4+(-4*j-6)*i^3+(6*j^2+6*j+11)*i^2+(-4*j^3+6*j^2-10*j-6)*i+j^4-6*j^3+11*j^2-6*j)*n^2+(-i^4+(16*j+6)*i^3+(-30*j^2-18*j-11)*i^2+(16*j^3-18*j^2+34*j+6)*i-j^4+6*j^3-11*j^2+6*j)*n-12*i*j*(i^2-3*i*j+j^2+1))*factorial(i+j-4)*(n-2)*factorial(2*n-i-4-j)*factorial(n)^2)/(factorial(n-j)*factorial(n-i)*factorial(2*n-3)*factorial(j)*factorial(i)); 
		aux3(i+1,j+1) = f_aux3(i,j,deg);
				
		% symmetric
		aux1(j+1,i+1) = aux1(i+1,j+1);
		aux2(j+1,i+1) = aux2(i+1,j+1);
		aux3(j+1,i+1) = aux3(i+1,j+1);
	end
	end
	
	% Evaluation of each ijop entry
	stiff = zeros(deg+1,deg+1,deg+1,deg+1);
	for i = 1:deg+1
	for j = 1:deg+1
		for o = 1:deg+1
		for p = 1:deg+1
			%stiff(i,j,o,p) = aux1(i,o)*aux2(j,p)+2*aux2(i,o)*aux2(j,p)+aux3(i,o)*aux1(j,p); 
			stiff(i,j,o,p) = aux1(i,o)*aux3(j,p)+2*aux2(i,o)*aux2(j,p)+aux3(i,o)*aux1(j,p);
		end
		end
	end
	end
	% Assemble the local stiffness matrix for one patch in a matrix.
	% Control points are numbered column after column : 
	% 		b00 b01 b02 ... b0n b10 b11 ...
	% Wouldn't it be, in the 1:M direction then the 1:N direction ?
	% Actually like matlab does ;-).
	% 		b00 b10 b20 ... bn0 , b01 b11 b21 ... bn1, b02, b12,...
	stiffLoc = reshape(stiff,(deg+1)^2,(deg+1)^2);
	
	% Assemble total stiffness matrix.
	% Patches are numbered column after column, so along 1:M direction, then along 1:N direction.
	stiff = zeros(ncp);
	for i = 1:M*N
		stiff((i-1)*(deg+1)^2+1:i*(deg+1)^2,(i-1)*(deg+1)^2+1:i*(deg+1)^2) = stiffLoc;
	end
	
	
	
	
	
	% ==================================================================  
	% Matrix B (dataMat)
	% ------------------
	% Assemble data term mass matrix from | Bx - d |^2
	% ==================================================================
	
	Bernstein = @(n,i,t) nchoosek(n,i)*t^i*(1-t)^(n-i);
	
	dataMat = zeros(ndp,ncp);
	for i = 1:ndp
		coeffs = zeros(deg+1);
		
		% find the closest patch of the data point.
		patchRow = max(1,min(M,ceil(dataCoords(i,1))));
		patchCol = max(1,min(N,ceil(dataCoords(i,2))));
		
		% coords of data point for patch m,n, between 0 and 1.
		locCoord = dataCoords(i,:) + 1 - [patchRow patchCol];
		
		for o = 1:deg+1
		for p = 1:deg+1
			coeffs(o,p) = Bernstein(deg,o-1,locCoord(1))*Bernstein(deg,p-1,locCoord(2));
		end
		end
		
		shift = (deg+1)^2*((patchCol-1)*M+patchRow-1)+1  :  (deg+1)^2*((patchCol-1)*M+patchRow);
		dataMat(i,shift) = coeffs(:);
	end
  
	
	
	% ==================================================================  
	% Matrix C.
	% ---------
	% Assemble the matrix of C1-continuity constraints.
	% ==================================================================
	
	% Assemble C1-continuity constraint matrices, for each interface.
	% Therefore, we construct four matrices:
	% 	- patchL : is the influence of a vertical interface on the
	% 			   patch on its left (continuity conditions on the right);
	% 	- patchR : is the influence of a vertical interface on the
	% 			   patch on its right (continuity conditions on the left);
	% 	- patchU : is the influence of a horizontal interface on the
	% 			   patch above (continuity conditions at the bottom);
	% 	- patchD : is the influence of a horizontal interface on the
	% 			   patch below (continuity conditions at the top);
	
	% Vertical interfaces
	patchL = zeros(2*(deg+1),(deg+1)^2);
	patchR = patchL; 					
	for i = 1:deg+1
		patchL(2*i-1,(deg-1)*(deg+1)+[i deg+1+i]) = [.5 -1];		% diff constraint
		patchL(2*i,  (deg-1)*(deg+1)+ i         ) =  .5;			% end of the diff constraint
		patchR(2*i,  [i deg+1+i]) = [-1 .5];						% diff constraint
		patchR(2*i-1,   deg+1+i)  =      .5;						% end of the diff constraint
	end
	
	% Horizontal interfaces
	patchU = zeros(2*(deg+1),(deg+1)^2);
	patchD = patchU;
	for i = 1:deg+1
		patchU(2*i-1,i*(deg+1)-[1 0]) = [.5 -1];					% diff constraint
		patchU(2*i,i*(deg+1)-1) = .5;								% interp constraint
		patchD(2*i-1,(i-1)*(deg+1)+2) = .5;							% interp constraint
		patchD(2*i,(i-1)*(deg+1)+[1 2]) = [-1 .5];					% diff constraint
	end
  
  
	% Assemble C1-continuity constraint matrix for all patches together
	% 	* N-1 vertical interfaces, times M patches,
	% 	* M-1 horizontal interfaces, times N patches, 
	% 	* 2*(deg+1) conditions by interfaces,
	% So there are 2*(deg+1)*(M(N-1) + N(M-1)) C1-conditions in total,
	% on the M*N*(deg+1)^2 control points.
	constraints = zeros(2*(deg+1)*(M*(N-1)+(M-1)*N),ncp);
	
	% M*(N-1) vertical constraints (patches L and R). 
	consInd = 0;
	for i = 1:M
	for j = 1:N-1
		% ranges
		consRange = consInd+1:consInd+2*(deg+1); 						% 2*(deg+1) constraints per interface
		pointsRangeL = (deg+1)^2*((j-1)*M+i-1)+1:(deg+1)^2*((j-1)*M+i); % points of the patch (i,j)
		pointsRangeR = (deg+1)^2*(j*M+i-1)+1:(deg+1)^2*(j*M+i);			% points of the patch (i,j+1)
		
		% Store the L and R constraints
		constraints(consRange,[pointsRangeL pointsRangeR]) = [patchL patchR];
		consInd = consInd + 2*(deg+1);
	end
	end
	
	% N*(M-1) horizontal constraints (patches U and D).
	for i = 1:M-1
	for j = 1:N
		% ranges
		consRange = consInd+1:consInd+2*(deg+1);						  % 2*(deg+1) constraints per interface
		pointsRangeU = (deg+1)^2*((j-1)*M+i-1)+1:(deg+1)^2*((j-1)*M+i);	  % points of the patch (i,j)
		pointsRangeD = (deg+1)^2*((j-1)*M+i)  +1:(deg+1)^2*((j-1)*M+i+1); % points of the patch (i+1,j)
		
		% Store the U and D constraints
		constraints(consRange,[pointsRangeU pointsRangeD]) = [patchU patchD];
		consInd = consInd+2*(deg+1);
	end
	end
	
	
	
	% ====================================================================  
	% Solve the linear system induced by minimizing the Euclidean 
	% quadratic fitting problem.
	% ====================================================================
  
	% prob : mat * x = ind * d
	%         [ A + lambda B'B     C' ] x = [ lambda B' ] d
	%         [       C            0  ]     [     0     ]
	mat = [ stiff  +  lambda*dataMat'*dataMat     constraints'; ...
				  constraints                   zeros(size(constraints,1))];

	ind = [lambda*dataMat'; ...
		 zeros(size(constraints,1),size(dataMat,1))];
			  
	% x = W d, so W = mat \ ind, normed.
	sol = mat \ ind;

	% Bezier control points now are first M*N*(deg+1)^2 entries of sol*data points
	sol = sol(1:ncp,:);
	sol(abs(sol)<thresh) = 0; 				 % set too small entries to zero.
	for i = 1:size(sol,1)
		sol(i,:) = sol(i,:) / sum(sol(i,:)); % Norm the solution such that each line is still an affine combination.
	end
  
      
	% ====================================================================
	% compute all control points - basically here is the weighted mean
	% b(i,j,m,n) = \sum_k alpha(i,j,m,n,k)
  % ====================================================================
	for m = 1:M	% all patches
	for n = 1:N
    for i = 1:deg+1	% each control point of the patch
    for j = 1:deg+1
      % the control point here under is a weighted mean of data
      % points mapped to the tangent space of the corner l.
      b{m,n,i,j} = 0;
      
      % Weighted sum over the dataPoints
      for k = 1:ndp
        weight = sol(((n-1)*M+(m-1))*(deg+1)^2+(j-1)*(deg+1)+i,k);
        if weight ~= 0
          b{m,n,i,j} = b{m,n,i,j} + weight * dataPoints{k}; % weighted mean
        end
      end
    end	% end each control point of the patch
    end
	end	% end all patches
	end


  % ====================================================================
  % Reshape the control points
  % ====================================================================
  controlPoints = cell(M*deg+1,N*deg+1);
  for m = 1:M
  for n = 1:N
    shiftX = (m-1)*deg;
    shiftY = (n-1)*deg;
    for i = 1:deg+1
    for j = 1:deg+1
      controlPoints{shiftX+i,shiftY+j} = b{m,n,i,j};
    end
    end
  end
  end
  
end











% function y = aux1(i,j,n).
% 	returns the value of the integral aux1(i,j,n) = int_[0,1] Bin(t)*Bjn(t) dt,
% 	where Bin and Bjn are ith and jth Bernstein polynomials of degree n.
function y = f_aux1(i,j,n)
	syms x;
	Bin = nchoosek(n,i)*(x^i)*((1-x)^(n-i));
	Bjn = nchoosek(n,j)*(x^j)*((1-x)^(n-j));
	
	y = double(int(Bin*Bjn,x,1,0));
end

% function y = aux2(i,j,n).
% 	returns the value of the integral aux1(i,j,n) = int_[0,1] dBin(t)*dBjn(t) dt,
% 	where dBin and dBjn are the derivatives of the ith and jth Bernstein 
% 	polynomials of degree n.
function y = f_aux2(i,j,n)
	syms x;
	dBin = diff(nchoosek(n,i)*(x^i)*((1-x)^(n-i)),x);
	dBjn = diff(nchoosek(n,j)*(x^j)*((1-x)^(n-j)),x);
	
	y = double(int(dBin*dBjn,x,1,0));
end

% function y = aux3(i,j,n).
% 	returns the value of the integral aux1(i,j,n) = int_[0,1] ddBin(t)*ddBjn(t) dt,
% 	where ddBin and ddBjn are the second derivatives of the ith and jth 
% 	Bernstein polynomials of degree n.
function y = f_aux3(i,j,n)
	syms x;
	ddBin = diff(diff(nchoosek(n,i)*(x^i)*((1-x)^(n-i)),x),x);
	ddBjn = diff(diff(nchoosek(n,j)*(x^j)*((1-x)^(n-j)),x),x);
	
	y = double(int(ddBin*ddBjn,x,1,0));
end

function pb = control_points_blended_2d(pb)
% Returns the control points of an approximating Bezier surface fitting
% a set of data points pb.dataPoints(i) at times pb.dataCoords(i) := (t1,t2).
%
% function pb = control_points_blended_2d(pb)
%    returns bp, the structure of the problem, updated with the control
% 	 points pb.controlPoints, the control points in their respective
%    tangent spaces (for blending purpose) pb.tangentPoints and 
%    the corner points of the patches pb.cornerPoints, which are also 
%    the root points of each tangent space.
%
%    The control points are obtained by a simple weighted mean done on
%    the tangent spaces of the corner points, by solving the Euclidean
%    minimization problem:
% 		min  ( 1/2 x^T A x + lambda/2 |Bx-d|^2 )
%        x
%		     s.t. Cx=0, 
% 	 where A is a stiffness matrix \int_[0,1]^2 < \D^2 B(t) / dt^2 > dt^2,
%          B is the Bezier matrix
%   	   d are the data points at time (t1,t2)
%          C is the matrix of C^1 global constraints
%          x are the unknown control points stacked in a vector patch per patch.
%
% Input: problem, the structure of the approximation problem, containing
% 			pb.dataPoints (the data points)
%   		pb.dataCoords (the coordinates (t1,t2) of each data point)
% 			pb.lambda     (the regularization parameter)
% 			pb.degree 	  (the degree of the curve)
% 		  pb.threshold  (the threshold on the weights of the final
% 	                  combination)
%       pb.patchM     (the number of patches in the first coordinate
%                     - rows, counted from top to bottom)
% 			pb.patchN     (the number of patches in the second coordinate
% 					          - columns, counted from left to right)
% 			pb.M	       (the manifold structure from Manopt)
%                   
% Output: problem, the structure of the approximation problem, edited
%         with two new fields named
%   		pb.tangentPoints (the (deg+1)^2 control points of each patch,
% 							  evaluated on the tangent spaces of the 
%                 four corners of the patch (m,n). They are
%                 stored in a (M,N,4,deg+1,deg+1) cell),
% 			pb.cornerPoints  (the (M+1,N+1) corners of the patches,
%                 directly given on the manifold).
%       pb.controlPoints (the (deg+1)^2 control points of each patch,
%                 stored in a (M+1,N+1)-cell and evaluated on the manifold).
%
%
% Original author: 
% 	Benedikt Wirth, Jan. 31, 2018.
% Contributors: 
%	Pierre-Yves Gousenbourger, Apr. 23, 2018
% Change log:
% 	Apr. 23, 2018 (PYG) - Documentation and integration of manopt.
% 	Apr. 26, 2018 (PYG) - Correction of indexes in the construction of C


 
% Important notes: 
%   - Bezier patches are counted from top left to bottom right, column 
%     by column. M is the number of rows and N the number of columns.
%   - Control points inside a Bezier patch are counted in the same way.
%   - The t1-coordinate points down (rows), the t2-coordinate points right (columns)

	% Data
	M = pb.patchM;
	N = pb.patchN;
	deg = pb.degree;
	lambda = pb.lambda;
	dataCoords = pb.dataCoords;
	dataPoints = pb.dataPoints;
	thresh = pb.threshold;
	man = pb.M;
	
	ndp = size(dataCoords,1); % number of data points
	ncp = M*N*(deg+1)^2; 	  % number of control points
	
	
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
	% compute corner control points of each patch - with the weighted mean
	% B(m,n) = exp( \sum_k alpha(...,k) log_{ref} (d_k)), on M,
	% where alpha is actually the "sol" vector.
	% ====================================================================
  
	%global_variables;
	g = @(x,y,t) man.exp(x,man.log(x,y),t);

	% Preallocation for control points
	B = cell(M+1,N+1); 				% B(m,n) is the top left control point of patch (m,n)
	b = cell(M,N,4,deg+1,deg+1); 	% b(m,n,l,i,j) is the 
									% 		(i,j)-control point in the 
									%  		(m,n)-Bezier-patch, expressed as
									% the log wrt the 
									% 		lth corner of the patch
									% (counting column-wise)
  cp = cell(M*deg+1,N*deg+1);  % cp contains all the control points on 
                  % the manifold. The reference point is the closest to
                  % the cp.
	logsD = cell(ndp,ndp); 
									% logs(h,k) is the logarithm of the 
									% 		hth data point wrt the 
									% 		kth data point
	logs = cell(ndp,M+1,N+1); 
									% logs(k,m,n) is the logarithm of the 
									% 		kth data point wrt the top left 
									% control point of the (m,n) patch
  
  
	for m = 1:M+1 % each patch
	for n = 1:N+1
	  
		% index of the control point
		controlPtIdx = ((n-1)*M+m-1)*(deg+1)^2+1;
		if m == M+1
			controlPtIdx = controlPtIdx - (deg+1)^2 + deg;
		end
		if n == N+1
			controlPtIdx = controlPtIdx - M*(deg+1)^2 + deg*(deg+1);
		end

		% The reference point to compute the control point should be the
		% data point with the largest weight.
		[~,refDataPtIdx] = max(abs(sol(controlPtIdx,:))); 

		% Weighted mean
		% B(m,n) = exp( \sum_k alpha(...,k) log_{ref} (d_k)), on M
		vec = 0;
		for k = 1:ndp
			if k ~= refDataPtIdx && sol(controlPtIdx,k) ~= 0
				if isempty(logsD{k,refDataPtIdx})
					%logsD{k,refDataPtIdx} = geo_log(dataPoints{refDataPtIdx},dataPoints{k});
					logsD{k,refDataPtIdx} = man.log(dataPoints{refDataPtIdx},dataPoints{k});
				end
				vec = vec + sol(controlPtIdx,k) * logsD{k,refDataPtIdx};
			end
		end
		%B{m,n} = geo_exp( dataPoints{refDataPtIdx}, vec );
		B{m,n} = man.exp( dataPoints{refDataPtIdx}, vec );
	  
	end	% end each patch
	end
      
  
      
	% ====================================================================
	% compute all control points - basically here is the weighted mean
	% b(i,j,m,n) = \sum_k alpha(i,j,m,n,k) log_{ref} (d_k), on T_{ref} M
	% The reference points are the control points on the corners of the
	% patch, that were computed in the previous section.
	% ====================================================================
	for m = 1:M	% all patches
	for n = 1:N
	  
		for l1 = 1:2   % all corner tangent spaces
		for l2 = 1:2
		  
			l = (l2-1)*2+l1;  	% numbering the corners from 1 to 4 as follows
                          % 1---3
                          % :   :
                          % 2---4
      
      corner = B{m+l1-1,n+l2-1};
			for i = 1:deg+1	% each control point of the patch
			for j = 1:deg+1
				% the control point here under is a weighted mean of data
				% points mapped to the tangent space of the corner l.
				b{m,n,l,i,j} = 0;
			  
				for k = 1:ndp	% each data point
			  
					% the weight of the data point is obtained within the sol matrix
					coeff = sol(((n-1)*M+(m-1))*(deg+1)^2+(j-1)*(deg+1)+i,k);
					%coeff = sol(((n-1)*M+N-1)*(deg+1)^2+(j-1)*(deg+1)+i,k);
				
					if coeff ~= 0
						if isempty(logs{k,m+l1-1,n+l2-1}) % map dataPoints{k} to the TP of the corner
							%logs{k,m+l1-1,n+l2-1} = geo_log(B{m+l1-1,n+l2-1},dataPoints{k});
							%logs{k,m+l1-1,n+l2-1} = man.log(corner,dataPoints{k});
              logs{k,m+l1-1,n+l2-1} = man.log(B{m+l1-1,n+l2-1},dataPoints{k});
						end
						b{m,n,l,i,j} = b{m,n,l,i,j} + coeff * logs{k,m+l1-1,n+l2-1}; % weighted mean
					end

				end	% end each data point
			  
			end	% end each control point of the patch
			end
      
      shiftx = (m-1)*deg+1+(l1-1)*2;
      shifty = (n-1)*deg+1+(l2-1)*2;
      
      posi = (l1-1)*2 + 1;
      posj = (l2-1)*2 + 1;
      
      cp{shiftx  ,shifty  } = man.exp( corner , b{m,n,l,posi  ,posj  });
      cp{shiftx+1,shifty  } = man.exp( corner , b{m,n,l,posi+1,posj  });
      cp{shiftx  ,shifty+1} = man.exp( corner , b{m,n,l,posi  ,posj+1});
      cp{shiftx+1,shifty+1} = man.exp( corner , b{m,n,l,posi+1,posj+1});
      
		end	% end all corners
		end
    
	end	% end all patches
	end

	pb.cornerPoints  = B;
	pb.tangentPoints = b;
  pb.controlPoints = cp;

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

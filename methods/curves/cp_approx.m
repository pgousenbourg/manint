function [cp, Tcp] = cp_approx(problem,lambda)
% Returns the control points of an approximating Bezier curve fitting
% the data of the problem.
%
% function cp = cp_approx(problem,lambda)
%    returns cp, a matrix of size [m,n,p] containing the control points
%    of the Bezier spline B approximating the data points problem.data.
%    lambda is the coefficient of fitting to data points.
%        min ||ddot(B)||^2 + lambda * sum || p_i - data_i ||^2
% 
% Input: problem, the structure of the approximation problem.
%        lambda, a scalar coefficient specifying the fitting to data
%                   
% Output: cp, the control points, stored in a [m,n,p] matrix where
%         [m,n] is the size of an element of the manifold problem.M.name
%         and p is the index of the control point of the Bezier curve.
%
% Original author: 
%   Pierre-Yves Gousenbourger, Oct. 13, 2016
% Contributors: 
%
% Change log:
% 	Oct. 13, 2016 (PYG):
%      First version
% 	Dec. 14, 2017 (PYG):
% 	   Modification for cubic fitting
% 	Jan. 23, 2018 (PYG) :
% 	   Modification of the calculus of C1 conditions (now done on the
% 	   tangent space). This modification (i) implies the storage of
% 	   the representation of points in the tangent space and their
%  	   root point ; (ii) overcomes the non-fitting curse observed for
%      periodic manifolds.
% 	Fev. 08, 2018 (PYG) :
% 	   Modification to have access to the representation of all control
% 	   points the tangent space of each data point.
% 	Fev. 22, 2018 (PYG) :
% 	   Modification of Tcp to be m x n x p x d, where d is the number of
% 	   data points. Tcp is sparse and stores the cp in the tangent
%      space of the data point i = 1...d.

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

	if ~isfield(problem,'reconstruction')
		type = 'LR';
	else
		type = problem.reconstruction;
	end

	M = problem.M;
    data = problem.data;
    [l,c,n] = size(data);
    
    
    % Prepare control points
	
	Nx   = 2*n;				 		% Number of unknowns
	Ncp  = 3*n-2;				 	% Number of control points at the end
	Tx  = zeros(l,c,Nx,n); 			% Vector of unknowns in tangent 
									%    spaces of all data points

    

    % Compute A and C in A*x = C*data
    
	A = sparse(Nx,Nx);
	C = sparse(Nx,n);
    
    switch n    
    case 2 ; cp = data ; Tcp = zeros(l,c,2,2); % 2 data points : cp = data.
	otherwise				
	% Case with 3 data points or more :
	% 	The A matrix and the C matrix can be computed by symetry as
	% 	the Bezier curve is reversible.
	%           			A*x = C*data
	% 	Note that, in the case with 3 data points, the loop is
	% 	never entered.
	
	
		% First segment : part due to the acceleration term
		A(1,1:4) = 3*[8 -12 2 2];
		A(2,1:4) = 3*[-12 24 -12 0];
		A(3,1:4) = 3*[2 -12 14 -4];
		A(4,1:4) = 3*[2 0 -4 2];
		
		% Inner segments : part due to the acceleration term
		for i = 2:n-2
			range = 2*i-1:2*i+2;
			A(2*i-1,range) = A(2*i-1,range) + 3*[2 -4 1 1];
			A(2*i  ,range) = A(2*i  ,range) + 3*[-4 14 -11 1];
			A(2*i+1,range) = A(2*i+1,range) + 3*[1 -11 14 -4];
			A(2*i+2,range) = A(2*i+2,range) + 3*[1 1 -4 2];
		end
		
		% Last segment : part due to the acceleration term
		A(end-3,end-3:end) = A(end-3,end-3:end) + 3*[2 -4 0 2];
		A(end-2,end-3:end) = A(end-2,end-3:end) + 3*[-4 14 -12 2];
		A(end-1,end-3:end) = A(end-1,end-3:end) + 3*[0 -12 24 -12];
		A(end  ,end-3:end) = A(end  ,end-3:end) + 3*[2 2 -12 8];
		
		% Part due to the fitting term
		A(1,1) = A(1,1) + 2*lambda;
		for i = 2:n-1
			A(2*i-1,2*i-1:2*i) = A(2*i-1, 2*i-1:2*i) + [1/2 1/2]*lambda;
			A(2*i  ,2*i-1:2*i) = A(2*i  , 2*i-1:2*i) + [1/2 1/2]*lambda;
		end
		A(end,end)           = A(end,end) + 2*lambda;
		
		
		
		% Computation of the C term
		C(1,1) = 2*lambda;
		for i = 2:n-1
			C(2*i-1,i) = lambda;
			C(2*i  ,i) = lambda;
		end
		C(end,end)     = 2*lambda;	
	
	
		% Matrix of weights
		D = inv(A)*C;
	

		% Compute the control points.
		% The C^1 conditions are matched on the tangent space T_di M as 
		% p_i = 0.5 * (b_i^+ + b_i^-).
		
		for i = 1:n
			
			% Projection of data on the tangent space of d
			d = data(:,:,i);
			Tdata = zeros(size(data));
			for k = 1:n
				Tdata(:,:,k) = M.log(d,data(:,:,k));
			end
			
			
			% Range of cp to compute for each method
			switch type
			case {'bezier','one','semi'};	range = 2*i-1:2*i;
			case 'LR';		range = max(1,2*i-3):min(Nx,2*i+2);
			case 'full';	range = 1:Nx;
			otherwise; error('Unknown type of generation.');
			end
					
			% Tx computation
			for j = range
				Tx(:,:,j,i) = cpGenW(D,j,Tdata);
			end
			
		end
		
		
		% C1 constraints and storage of points on T_di M.
		% Tcp is a matrix of size [n,m,Ncp,Ndata] where
		% 	[n,m] is the size of the embedding space ;
		% 	Ncp is the number of control points ;
		% 	Ndata is the number of data points.
		
		Tcp 	= zeros(l,c,Ncp,n);
		
		% Store Tx to Tcp
		Tcp(:,:,sort([1,  2:3:end-2,  3:3:end-1,  end]),:) = Tx;
		% C1 conditions
		Tcp(:,:,4:3:end-3,:) = 0.5*(Tcp(:,:,3:3:end-4,:) + Tcp(:,:,5:3:end-2,:));
		
		
		% Cleaning in LR
		if strcmp(type,'LR')
			for i = 2:n-1
				Tcp(:,:,3*(i-1)+2,i-1) = zeros(l,c);
				Tcp(:,:,3*(i-1),i+1)   = zeros(l,c);
			end
		end
		
		
		
		% Storage of the actual control points on M
		cp = zeros(l,c,Ncp);
		
		cp(:,:,1) = M.exp(data(:,:,1),Tcp(:,:,1,1));
		cp(:,:,2) = M.exp(data(:,:,1),Tcp(:,:,2,1));
		for i = 2:n-1
			idx = 3*(i-1)+1;
			cp(:,:,idx-1) = M.exp(data(:,:,i),Tcp(:,:,idx-1,i));
			cp(:,:,idx)   = M.exp(data(:,:,i),Tcp(:,:,idx  ,i));
			cp(:,:,idx+1) = M.exp(data(:,:,i),Tcp(:,:,idx+1,i));
		end
		cp(:,:,end-1) = M.exp(data(:,:,end),Tcp(:,:,end-1,end));
		cp(:,:,end)   = M.exp(data(:,:,end),Tcp(:,:,end,end));
	end
end





% Auxiliary function to compute the weighted sums
function y = cpGenW(D,j,Tdata)
	[a,b,~,~] = size(Tdata);
	% weights
	w = reshape(full(D(j,:)),1,1,length(D(j,:)));
	w = repmat(w,[a,b,1]);

	% weighted sum
	y = sum(w.*Tdata,3);
end

function [cp, Tcp] = cp_interp(problem)
% Returns the control points of an interpolating Bezier curve fitting
% the data of the problem.
%
% function cp = cp_interp(problem)
%    returns cp, a matrix of size [m,n,p] containing the control points
%    of the Bezier spline B interpolating the data points problem.data,
% 	 such that
%        min ||ddot(B)||^2 
% 
% Input: problem, the structure of the approximation problem.
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
% 	Mar. 22, 2018 (PYG):
%      First version

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
	
	Nx 	 = n;						% Number of unknowns
	Ncp  = 3*n-2;				 	% Number of control points at the end
	Tx  = zeros(l,c,Nx,n); 			% Vector of unknowns in tangent 
									% spaces of all data points
	Tcp 	= zeros(l,c,Ncp,n);		% Tcp is a matrix of size [n,m,Ncp,Ndata]
									% that stores the control points
									% in the tangent space T_di M
		

    

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
	
	
		% First segment
		A(1,1:2) = [6 -3];
		A(2,1:2) = [-3 6];
		
		C(1,1) = 3;
		C(2,2) = 3;
		
		% Other segments
		for i = 2:n-1
			A(i  ,i:i+1) = A(i  ,i:i+1) + [6 3];
			A(i+1,i:i+1) = A(i+1,i:i+1) + [3 6];
			
			C(i  ,i)     = C(i  ,i    ) + 9;
			C(i+1,i:i+1) = C(i+1,i:i+1) + [6 3];
		end
			
	
		% Matrix of weights
		D = inv(A)*C;
	

		% Compute the control points.
		% The C^1 conditions are matched on the tangent space T_di M as 
		% b_i^+ = 2*d_i - b_i^-
		
		for i = 1:n
			
			% Projection of data on the tangent space of d
			d = data(:,:,i);
			Tdata = zeros(size(data));
			for k = 1:n
				Tdata(:,:,k) = M.log(d,data(:,:,k));
			end
			
			% Range of cp to compute for each method
			switch type
			case {'bezier','one','semi'};	range = max(1,i-1):min(Nx,i+1);
			case 'LR';		range = max(1,i-1):min(Nx,i+1);
			case 'full';	range = 1:Nx;
			otherwise; error('Unknown type of generation.');
			end
					
			% Tx computation
			for j = range
				Tx(:,:,j,i) = cpGenW(D,j,Tdata);
			end
			
			% Store the data points
			Tcp(:,:,1:3:end,i) = Tdata;
		end
		
		
		% C1 constraints and storage of points on T_di M.
		% Store Tx to Tcp
		Tcp(:,:,[2,3:3:end-1],:) = Tx;
		% C1 conditions
		Tcp(:,:,5:3:end-2,:) = 2*Tcp(:,:,4:3:end-3,:) - Tcp(:,:,3:3:end-4,:);
		
		
		%% Cleaning in LR
		%if strcmp(type,'LR')
			%for i = 2:n-1
				%Tcp(:,:,3*(i-1)+2,i-1) = zeros(l,c);
				%Tcp(:,:,3*(i-1),i+1)   = zeros(l,c);
			%end
		%end
		
		
		
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

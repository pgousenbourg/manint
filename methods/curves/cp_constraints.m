function cp = cp_constraints(x,Ncp,M)
% Returns the control points of an approximating Bezier curve fitting
% the data of the problem.
%
% function cp = cp_constraints(x,Ncp)
%    Returns the control points cp by applying the constraints of the
%    Bezier problem to x.
% 
% Input: x - the solution of the optimization
%        Ncp - the number of control points
% 		 M - the manifold
%                   
% Output: cp, the control points, stored in a [m,n,p] matrix where
%         [m,n] is the size of an element of the manifold problem.M.name
%         and p is the index of the control point of the Bezier curve.
%
% Original author: 
%   Pierre-Yves Gousenbourger, Jan. 23, 2017
% Contributors: 
%
% Change log:
% 	Jan. 23, 2017 (PYG):
%      First version
% 	Dec. 14, 2017 (PYG):
% 	   Modification for cubic fitting.

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

	Ndata = (Ncp + 2)/3;
	cp = zeros(size(x,1),size(x,2),Ncp); 
	
	% bm and bp
	cp(:,:,2:3:end-2) = x(:,:,2:2:end-2);   % bp
	cp(:,:,3:3:end-1) = x(:,:,3:2:end-1);   % bm
	
	% First and last control points
	cp(:,:,1) 	= x(:,:,1);					% p0
	cp(:,:,end) = x(:,:,end);				% pn
	
	for i = 2:Ndata-1
		cp(:,:,3*(i-1)+1) = M.exp(cp(:,:,3*(i-1)),0.5*M.log(cp(:,:,3*(i-1)), cp(:,:,3*(i-1)+2)));   % p_i
	end
end

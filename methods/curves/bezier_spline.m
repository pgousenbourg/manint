function y = bezier_spline(problem,t,display)
% Returns the value of the Bezier spline controlled by control points cp
% at time t.
%
% function y = bezier_spline(problem,t,method)
%    returns y, a matrix of size [m,n,p] containing the value of the
%    hybrid Bezier spline B at times t. B is composed of second order 
%    Bezier functions as first and last segments, and cubic Bezier 
%    functions as intermediate segments. t may be a vector of values
%    with 0 <= t <= (size(cp) + 2)/3
% 
% Input: problem, the structure of the problem with 
% 			control : control points on M
% 			Tcp 	: control points on the tangent spaces of di
% 			data 	: the data points
% 		 t, the time at which the curve is evaluated
% 		 display, if a waitbar must be displayed
%                  
% Output: y, the value of the spline at t, stored in a [m x n x p].
%
% Original author: 
%   Pierre-Yves Gousenbourger, Oct. 19, 2016
% Contributors: 
%
% Change log:
% 	Oct. 19, 2016 (PYG):
%      First version
% 	Dec. 14, 2017 (PYG):
% 	   Modification for cubic splines

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
		display = 'false';
	end
	wb = strcmp(display,'true');
		
	% Data
	cp 		= problem.control;				% computed control points
	M  		= problem.M;					% manifold
	[l,c,m] = size(cp);                     % (l,c), size of the space
											% m, number of cp
	n 	  	= (m + 2)/3;					% number of data points
	p 		= length(t);                    % number of computed points
	f = @(x,ti) decasteljau(M,x,ti);		% method for reconstruction
	
	
	
	% Waitbar because it is fun...
	if wb; h = waitbar(0,'Reconstruction. Please wait...');  total = p; growth = 0; end
	
	
	% Preallocation
	y 		= zeros(l,c,p);
		


	% Compute the Bezier spline.
	% In the case of only two control points, only a geodesic is needed.
	% In the other cases, we apply f.
	switch m
	case 2
		a = cp(:,:,1);
		b = cp(:,:,2);
		for i = 1:p
			y(:,:,i) = M.exp(a,M.log(a,b),t(i));
			if wb; growth = growth + 1;	 waitbar(growth/total,h,sprintf('Reconstruction. Complete: %.f %%.',growth/total*100)); end
		end

	otherwise	
		for i = 1:p
			ti   = t(i);
			sgmt = floor(ti)+1;
			
			if sgmt == n ; y(:,:,i) = cp(:,:,3*n-2); % last point
			elseif sgmt == ti+1 ; y(:,:,i) = cp(:,:,3*(sgmt-1)+1); % first or inner points
			else
				range = [3*(sgmt-1)+1:3*sgmt+1];
				ti = ti - sgmt + 1;
				
				y(:,:,i) = f(cp(:,:,range),ti);
			end
			if wb; growth = growth + 1;	 waitbar(growth/total,h,sprintf('Reconstruction. Complete: %.f %%.',growth/total*100)); end
		end
	end
	
	% End of waitbar
	if wb; close(h); end
end

function y = blended_spline(problem,t,display)
% Returns the value of the Bezier spline controlled by control points cp
% at time t.
%
% function y = blended_spline(problem,t)
%    returns y, a matrix of size [m,n,p] containing the value of the
%    cubic Bezier spline B at times t.
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
% 	Jan. 24, 2018 (PYG):
% 	   Modification of the reconstruction process. The first step is now
% 	   an Exp mapping of the straight line linking b_i^+, p_i and b_i^-
% 	   in the tangent space of d_i. Therefore, the curve is no more a
% 	   Bezier curve, but a Bezier-like curve.
% 	Fev. 27, 2018 (PYG) :
% 		Modification of the reconstruction process. On the i^th segment,
% 		the De Casteljau algorithm is fully performed in the tangent
% 		space of p_i and p_{i+1}. Then, a weighted mean is done on the
% 		two obtained points.

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
		display = 'true';
	end
	wb = strcmp(display,'true');
	
	% parameters
	
	d	 	  = problem.data;					% data points on M
	Tcp 	  = problem.Tcp;					% cp on T_di M
	M  		  = problem.M;						% manifold
	[l,c,m,n] = size(Tcp);                     	% (l,c), size of the space
												% m, number of cp
												% n, number of data points
	p 		  = length(t);                    	% size of the discretization
	f 		  = @(x,ti) euclideanBezier(x,ti);	% reconstruction method
	
	
	% Waitbar because it is fun...
	if wb; h = waitbar(0,'Reconstruction. Please wait...'); total = p; growth = 0; end;
	


	% Reconstruction of the curve
	y 		= zeros(l,c,p);
	
	
	switch m
	case 2		% 2 cp: geodesic
		a = d(:,:,1);  b = d(:,:,2);
		
		for i = 1:p 
			y(:,:,i) = M.exp(a,M.log(a,b),t(i));
			
			% Waitbar
			if wb
				growth = growth + 1; 
				waitbar(growth/total,h,sprintf('Reconstruction. Complete: %.f %%.',growth/total*100));
			end
		end



	% more than 2 cp : three steps
	% 	1- Compute the Bezier spline on T pi M and T p(i+1) M.
	% 	2- Map the result to M.
	% 	3- Weighted mean of the two points.
	otherwise		
		x = zeros(l,c,2);
		
		for i = 1:p
			ti   = t(i);
			sgmt = floor(ti)+1;
			
			if sgmt == n ; y(:,:,i) = M.exp(d(:,:,n),Tcp(:,:,3*n-2,n)); % last point
			elseif sgmt == ti+1 ; y(:,:,i) = M.exp(d(:,:,sgmt),Tcp(:,:,3*(sgmt-1)+1,sgmt)); % first or inner points
			else
				range = [3*(sgmt-1)+1:3*sgmt+1];
				
				% De casteljau on the two tangent spaces
				ti = ti - sgmt + 1;
				x(:,:,1) = M.exp(d(:,:,sgmt)  ,f(Tcp(:,:,range,sgmt  ),ti));
				x(:,:,2) = M.exp(d(:,:,sgmt+1),f(Tcp(:,:,range,sgmt+1),ti));
				
				w = weighting(ti);
				% Weighted sum
				y(:,:,i) = weighted_sum([1-w w],x,M);
			end
			if wb; growth = growth + 1;	 waitbar(growth/total,h,sprintf('Reconstruction. Complete: %.f %%.',growth/total*100)); end
		end
	end
	
	% End of waitbar
	if wb; close(h); end
end

% Auxiliary function to compute the weights
function s = weighting(ti)
	s = 3*ti^2 - 2*ti^3;
end


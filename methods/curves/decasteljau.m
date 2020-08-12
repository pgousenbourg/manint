function y = decasteljau(M,b,t)
% Returns the value of a Bezier function driven by control points b at
% time t.
%
% function y = decasteljau(M,b,t)
%    returns y, the value of the manifold valued Bezier function on M,
%    controlled by the control points b of M, and at time t (scalar). 
%    y is a vector of the dimension of the manifold.
% 
% Input: M, the manifold structure (from manopt framework) containing
%         the elements of differential geometry.
%        b, the control points of the Bezier function (evaluated in M).
%        t, the scalar time at which the spline is evaluated.
%                  
% Output: y, the value of the spline at t.
%
% Original author: 
%   Pierre-Yves Gousenbourger, Oct. 19, 2016
% Contributors: 
%
% Change log:
% 	Oct. 19, 2016 (PYG):
%      Firts version

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
    
    assert(t <= 1 & t >= 0,'t must be between [0,1]');
    
    % tool
    geo_map = @(x,y,t) M.exp(x,t*M.log(x,y));
    
    % For each step in the de casteljau process
    y = b;
	for i = 1:size(b,3)-1
		for j = 1:size(y,3)-1
			y(:,:,j) = geo_map(y(:,:,j),y(:,:,j+1),t);
		end
	end
	
	y = y(:,:,1);
end

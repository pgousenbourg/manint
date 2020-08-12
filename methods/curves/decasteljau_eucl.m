function y = decasteljau_eucl(b,t)
% Returns the value of a Euclidean Bezier function driven by Euclidean
% control points b at time t.
%
% function y = decasteljau_eucl(b,t)
%    returns y, the value of the Euclidean Bezier function,
%    controlled by the control points b, and at time t (scalar). 
%    y is a vector of the dimension of the space of b.
% 
% Input: b, the control points of the Bezier function.
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
%      	First version
% 	Feb. 08, 2018 (PYG):
% 	   	Modification to only euclidean

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
    
    % weighted mean
    f = @(x,y,t) (1-t)*x + t*y ;
    
    % For each step in the de casteljau process
    y = b;
	for i = 1:size(b,3)-1
		y = f(y(:,:,1:end-1),y(:,:,2:end),t);
	end
	
	y = y(:,:,1);
end

%% DE_CASTELJAU_EUCLIDEAN
%       Computes the Bezier surface in euclidean space based on the Bezier points b,
%       evaluated at time (t1,t2) in [0,1]. The degree (d) of the Bezier
% 		surface is given in entry.
% 
% Input: b       : the control points (cell d+1 x d+1).
%        t1,t2   : the evaluated time on the square [0,1]x[0,1].
% 		 d 		 : the degree of the curve.
%
% Output: the Bezier curve, in a matrix [dim1 x dim2]

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

function y = de_casteljau_euclidean(b,d,t1,t2)
  assert(size(b,1) <= d+1,'  The degree of the curve is too high for the number of control points provided');

	x = b;
	
	% De Casteljau
	for k = 1:d
		for i = 1:d-k+1
			for j = 1:d-k+1
				% Conrol Points for the average
				p = x(i:i+1,j:j+1);
				% Compute the average
                x{i,j} = (1-t1)*(1-t2)*p{1,1} + t1*(1-t2)*p{2,1} + (1-t1)*t2*p{1,2} + t1*t2*p{2,2};
            end
        end
        % Delete last entries
		x = x(1:end-1,1:end-1);
	end
	
	% Verify that size x = 1
	assert(size(x,1) == 1 && size(x,2) == 1);
	
	y = squeeze(x{:});
end

function w = blendSurfWeights(X,Y,rootCoords)
% Blending weights used in the blending technique of surfaces.
%
% function w = blendSurfWeights(X,Y)
%    returns w, a 4-cell containing the four weights used in the
%    blending of surfaces at times (X,Y) in [0,1]x[0,1]. By default, the
%    positions of the linearization points are supposed to be at the 
%    corners of this unitary patch.
%    Each weights entry is a matrix of the same size as x and y.
%    The weights are stored in the following (graphical) manner:
%        
%         01  ---  11       ->     w{3}  ---  w{4}
%          |       |                |          |
%          |       |                |          |
%         00  ---  10              w{1}  ---  w{2}
%
% function w = blendSurfWeights(X,Y,rootCoords)
%    returns w, a 4-cell containing the four weights used in the
%    blending of surfaces at time (X,Y) in the patch defined by the 
%    rootCoords.
%    When rootCoords is not specified, it falls down naturally to
%    the unitary patch [0,1]x[0,1].
%
% If X,Y are vectors with only one entry, then w is returned as a
% 4-vector.
%
% EXAMPLE: 
%    >> X = 0; Y = 0;
%    >> rootCoords = [0 0; 1 0; 0 1; 1 1];
%    >> blendSurfWeights(X,Y,rootCoords);
%    
%      ans = 
%         
%          1  0  0  0
%
% EXAMPLE 2:
%    >> X = [0 0 1 1];
%    >> Y = [0 1 0 1];
%    >> blendSurfWeights(X,Y)
%    
%      returns a 4-cell of [1x4]-matrices.
%
%
% inputs:  X,Y are coordinates of the same size ([pxq]-matrix);
%          rootCoords is a [4x2]-matrix of coordinates of the four 
%            corners of the patch.
%
% outputs: w is a [4]-cell of weights, where each entry is a 
%          [pxq]-matrix, except if p=q=1 (in that case, w is a
%          4-vector.
%
% Original author: 
% 	Pierre-Yves Gousenbourger, Dec. 18, 2019.
% Contributors: 
%
% Change log:
% 	Dec. 18, 2019 (PYG) - First version.
%   Jan. 06, 2020 (PYG) - Modification of the logic to stage better the
%                         with the rest of the project. Now weights are
%                         stored in 4-cells.

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

  % parameters  
  sX = size(X);
  lX = prod(sX);
  
  % create rootCoords if necessary
  if nargin == 2
    rootCoords = [0 0; 1 0; 0 1; 1 1];
  end

  % defense
  % x and y must be the same size
  defense = prod(size(X) == size(Y)) == 1;
  assert(defense,'surfaceFitting:dimCheck','X and Y must be of same size')
  
  % on input formatting: rootCoords must be of size 4x2
  assert(size(rootCoords,2) == 2 && size(rootCoords,1) == 4,'surfaceFitting:dimCheck','rootCoords must be a (4x2)-matrix');
  
  
  % computation of the weights
  for i = 1:2;
  for j = 1:2;
    k  = rootPos(i,j);
    t = abs([X(:),Y(:)] - rootCoords(k,:));
    
    % width of the windows in X and Y
    width = [0 0];
    width(1,1) = rootCoords(rootPos(2,j),1) - rootCoords(rootPos(1,j),1);
    width(1,2) = rootCoords(rootPos(i,2),2) - rootCoords(rootPos(i,1),2);
    width = repmat(width,[lX,1]);
    
    % defense
    assert(sum(sum((width == 0)))==0,'surfaceWeights:gridMisAlignment','The given rootCoords implied a 0-sized width.\nIs the grid turned by 90 degrees?');
    
    % weights for each entry reshaped at the size of X,Y 
    w{k} = reshape(prod(g(t./width),2),sX);
  end
  end
  
  % case where there is only one point
  if length(sX) == 2 && prod(sX) == 1
    w = cell2mat(w);
  end
  
end

% Auxiliary function of weights
function s = g(t)
	s = 2*t.^3 - 3*t.^2 + 1;
end

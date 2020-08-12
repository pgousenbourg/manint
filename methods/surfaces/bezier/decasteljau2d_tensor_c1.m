function y = decasteljau2d_tensor_c1(varargin)
% Computes a Bezier curve by double tensorization.
%
% function y = decasteljau2d_tensor_c1(M,b,t1,t2,patchX,patchY,maxPatchX,maxPatchY)
%   returns y, the value of the Riemannian Bezier curve on the manifold
%   M at t1 and t2. The curve is driven by the points b (cell).
%   This method is for C1-patching of Bezier surfaces only. It requires
%     patchX - the number of the current patch (in X)
%     patchY - the number of the current patch (in Y)
%     maxPatchX - the max number of patches in the surface (in X)
%     maxPatchY - the max number of patches in the surface (in Y)
%
% function y = decasteljau2d_tensor_c1(M,b,t1,t2,patchX,patchY,maxPatchX,maxPatchY,'direction',DIRECTION)
%   returns y, the value of the Riemannian Bezier curve on the manifold
%   M by applying first the tensorization in the direction DIRECTION.
%
% Original author: 
% 	Pierre-Yves Gousenbourger, Feb. 21, 2019.
% Contributors: 
%   ***
% Change log:
% 	Feb. 21, 2019 (PYG) - Original file.

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

  % create the parser
  ip = inputParser();
  addRequired(ip,'M');
  addRequired(ip,'b');
  addRequired(ip,'t1');
  addRequired(ip,'t2');
  addRequired(ip,'patchX');
  addRequired(ip,'patchY');
  addRequired(ip,'maxPatchX');
  addRequired(ip,'maxPatchY');
  addOptional(ip,'direction','x-direction'); % x-direction by default

  % parse the parser as a structure.
  parse(ip, varargin{:});
  vars = ip.Results;
  
  % variables
  M  = vars.M;
  b  = vars.b;
  t1 = vars.t1;
  t2 = vars.t2;
  patchX = vars.patchX;
  patchY = vars.patchY;
  maxPatchX = vars.maxPatchX;
  maxPatchY = vars.maxPatchY;
  direction = vars.direction; 

  % Check the direction
  if strcmp(direction,'y-direction'); 
    b = b';
    % switch times
    temp = t1; t1 = t2; t2 = temp;
    % switch patches
    temp = patchX; patchX = patchY; patchY = temp;
    % switch maxPatches
    temp = maxPatchX; maxPatchX = maxPatchY; maxPatchY = temp;
    
    clear temp;
  end
  
  % parameters
  deg = size(b,1);
  
  % fictional control points
  bb = cell(deg,1);
  for i=1:deg
    bb{i} = decasteljau1d_c1(M,b(:,i),t1,patchX,maxPatchX);
  end
  
  % solution
  y = decasteljau1d_c1(M,bb,t2,patchY,maxPatchY);

end


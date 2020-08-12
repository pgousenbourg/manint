function y = decasteljau2d_mean_c1(varargin)
% Computes a Bezier curve with a 2d-karcher-mean approach such that
% the Bezier surface is C1
%
% function y = decasteljau2d_mean_c1(M,b,t1,t2,patchX,patchY,maxPatchX,maxPatchY)
%   returns y, the value of the Riemannian Bezier curve on the manifold
%   M at t1 and t2. The curve is driven by the points b (cell).
%   This method is for C1-patching of Bezier surfaces only. It requires
%     patchX - the number of the current patch (in X)
%     patchY - the number of the current patch (in Y)
%     maxPatchX - the max number of patches in the surface (in X)
%     maxPatchY - the max number of patches in the surface (in Y)
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
  
  deg = size(b,1)-1;
  
  % Curve reconstruction
  w_x = weights_bern_c1(patchX,deg,t1,maxPatchX);
  w_y = weights_bern_c1(patchY,deg,t2,maxPatchY);
  w = w_x*w_y';
  
  % Manopt solution for the Karcher mean
  model.M = M;
  model.cost = @(x) karcher(M,x,b,w);
  model.grad = @(x) dkarcher(M,x,b,w);
  
  options.verbosity = 0;
  warning('off', 'manopt:getHessian:approx') 
  
  % Solution
  y = trustregions(model,[],options);
end


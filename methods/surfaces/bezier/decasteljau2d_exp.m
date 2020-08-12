function y = decasteljau2d_exp(varargin)
% Computes a Bezier curve by 2d-deCasteljau with exponential map.
% DISCLAIMER: THIS PROGRAM IS ONLY FOR PEDAGOGIGAL PURPOSES.
%
% function y = decasteljau2d_exp(M,b,t1,t2)
%   returns y, the value of the Riemannian Bezier curve on the manifold
%   M at t1 and t2. The curve is driven by the points b (cell).
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

  % parse the parser as a structure.
  parse(ip, varargin{:});
  vars = ip.Results;
  
  % variables
  M  = vars.M;
  b  = vars.b;
  t1 = vars.t1;
  t2 = vars.t2;
  
  w_x = [1-t1 ; t1];
  w_y = [1-t2 ; t2];
  w = w_x*w_y';
  
  deg = size(b);
  assert(deg(1)==deg(2),'The control points must be a square-cell.')
  deg = deg(1)-1;
  
  % loop
  for k = 1:deg
    for i = 1:deg-k+1;
    for j = 1:deg-k+1;
      % Lifting
      root = b{i,j};
      cp{1,1} = M.log(root,b{i,j});
      cp{2,1} = M.log(root,b{i+1,j});
      cp{1,2} = M.log(root,b{i,j+1});
      cp{2,2} = M.log(root,b{i+1,j+1});      
      
      % mean
      Tb = cp{1,1}*w(1,1) + cp{1,2}*w(1,2) + cp{2,1}*w(2,1) + cp{2,2}*w(2,2);
      
      % back-projection
      b{i,j} = M.exp(root,Tb);
    end
    end
  end
  
  % solution
  y = b{1,1};
end

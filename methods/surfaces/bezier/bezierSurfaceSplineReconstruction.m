function Z = bezierSurfaceSplineReconstruction(X,Y,controlPoints,rootCoords,degree)
% Returns the composite Euclidean Bezier surface of a given degree, driven 
% by the controlPoints and evaluated at X and Y.
%
% function Z = bezierSurfaceSplineReconstruction(X,Y,controlPoints,degree)
%    returns Z, the value of the Euclidean composite Bezier surface of 
%    degree DEGREE, driven by the control points CONTROLPOINTS and 
%    evaluated at (X,Y). The CONTROLPOINTS are stored in a (mcp*ncp)-cell
%    of matrices (m*n), where m*n is the size of the Euclidean space.
%    The values X,Y are matrices of same size [p*q]. The degree is a scalar.
%
% outputs: Z is a [p,q,m,n] matrix.
%
% Original author: 
% 	Pierre-Yves Gousenbourger, Mar. 30, 2015.
% Contributors: 
%	  Paul Striewski, Jul. 14, 2015.
%   Anna Helene Scholzen, Apr. 28, 2018.
% Change log:
%   Mar. 20, 2015 (PYG) - First version of Bezier surfaces reconstruction.
%   Apr. 28, 2018 (AHS) - Modification of the version to blended surfaces.
% 	Nov. 07, 2018 (PYG) - Integration of the Manopt framework to Bezier
%                         surfaces reconstruction.
%   May. 13, 2019 (PYG) - Merge of the versions of Apr. 28 and Nov. 07.
% 	May. 06, 2020 (PYG) - Alignment with the rest of the packaged code.

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
  if nargin < 5
    degree = 3;
  end
  patchM = (size(controlPoints,1)-1)/degree;
  patchN = (size(controlPoints,2)-1)/degree;
  
  % defense  
  assert(floor(patchM)==patchM && floor(patchN)==patchN,'The number of control points does not correspond to the degree of the curve');
  assert(ismatrix(controlPoints));
  
  % other parameters
  if isvector(X)
    sX = length(X);
  else
    sX = size(X);
  end
  X = X(:);
  Y = Y(:);
  assert(length(X) == length(Y),'surfaceFitting:dimCheck','X and Y must be of same size');
  sData  = size(controlPoints{1,1});
  
  
  % preallocation
  Z = zeros(prod(sX),sData(1),sData(2));
  
  % reparameterizing the [X,Y] w.r.t the dataCoords such that they belong
  % to the domain of the Bezier surface, where each patch is between [i,i+1]x[j,j+1].
  X = (X - min(rootCoords(:,1))).*patchM./(max(rootCoords(:,1)) - min(rootCoords(:,1)));
  Y = (Y - min(rootCoords(:,2))).*patchN./(max(rootCoords(:,2)) - min(rootCoords(:,2)));
  
  % Reconstruction for each value of X,Y
  for i = 1:length(X)
    % find in which patch we are
    patchX = max(1,ceil(X(i)));
		patchY = max(1,ceil(Y(i)));
    
    % reparametriez between 0 and 1
    t1 = X(i) - patchX + 1;
		t2 = Y(i) - patchY + 1;
    
    % get control points
    rangeX   = (patchX-1)*degree + 1 : (patchX)*degree + 1;
		rangeY   = (patchY-1)*degree + 1 : (patchY)*degree + 1;
    controls = controlPoints(rangeX,rangeY);
    
    % Bezier Surface Reconstruction
    Z(i,:,:) = bezierSurfaceEuclidean(controls,t1,t2);
  end
  
  % reshape the solution to the dimension the initial X and Y
  Z = reshape(Z,[sX,sData(1),sData(2)]);
  
end

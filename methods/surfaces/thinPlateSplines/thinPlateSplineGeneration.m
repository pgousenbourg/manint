function repr = thinPlateSplineGeneration(dataPoints,dataCoords,lambda)
% Returns the coefficients conducting the Euclidean thin Plate Spline.
%
% function repr = thinPlateSplineGeneration(dataPoints,dataCoords)
%    returns repr, a structure containing repr.a and repr.d,
%    the coefficients driving the Thin Plate Spline
%    interpolating the dataPoints{i} at dataCoords(i).
%
% function repr = thinPlateSplineGeneration(dataPoints,dataCoords,lambda)
%    returns repr, a structure containing repr.a and repr. d, 
%    the coefficients driving the Thin Plate Spline
%    fitting the dataPoints{i} at dataCoords(i) with a fitting parameter
%    lambda (the smaller, the closer).
%
% inputs:  dataPoints is a n-dimensional cell. Stores the data points
%            of size [p,q].
%          dataCoords is a [nx2]-matrix. Stores the coordinates at which
%            the data points must be fitted.
%          lambda     is a scalar (optional). The smaller, the closer to
%            interpolation the curve will be.
%
% outputs: repr. d is a [n,p,q]-dimensional matrix, collecting the values of 
%            the coefficients associated to the rbf.
%          repr. a is a [3,p,q]-dimensional matrix, with the coefficients 
%            associated to the affine part of the TPS.
%
% Original author: 
% 	Pierre-Yves Gousenbourger, Dec. 17, 2019.
% Contributors: 
%
% Change log:
% 	Dec. 17, 2019 (PYG) - First version.

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

  % Defense
  assert(iscell(dataPoints),'surfaceFitting:dataNotInCell','The data points must be stored in a cell');
  dataPoints = dataPoints(:);
  
  assert(size(dataCoords,2) == 2,'surfaceFitting:coordsNotIn2D','The dataCoords must be stored in a nx2 matrix');
  assert(size(dataCoords,1)==length(dataPoints),'surfaceFitting:dimCheck','There must be as many dataPoints as dataCoords');
  if nargin < 3
    lambda = 0;
  end
  
  % parameters
  ndp   = length(dataPoints);
  [m,n] = size(dataPoints{1});
  
  % Process dataPoints
  DP = reshape(cell2mat(dataPoints'),[m,n,ndp]);
  DP = reshape(permute(DP,[3,1,2]),[ndp,m*n]);
  
  % Construction of the matrix
  E = zeros(ndp,ndp);
  T = zeros(3,ndp);
  I = eye(ndp,ndp);
  
  for i = 1:ndp
    for j = 1:ndp
      r = norm(dataCoords(i,:) - dataCoords(j,:));
      E(i,j) = rbf(r);
    end
    T(:,i) = [1,dataCoords(i,:)]';
  end
  
  % Solve the problem
  x = [E+lambda*I T'; T zeros(3,3)]\[DP;zeros(3,m*n)];
  
  % extract a and d
  d = reshape(x(1:ndp,:),[ndp,m,n]);
  a = reshape(x(end-2:end,:),[3,m,n]);
  
  repr.a = a;
  repr.d = d;
end

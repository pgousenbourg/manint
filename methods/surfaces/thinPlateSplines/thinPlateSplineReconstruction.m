function tps = thinPlateSplineReconstruction(X,Y,repr,t)
% Returns the Euclidean Thin Plate Spline driven by a and d at times
% X and Y.
%
% function tps = thinPlateSplineReconstruction(X,Y,repr,t)
%    returns tps, the value of the Euclidean Thin Plate Spline driven 
%    by the coefficients repr.a and repd.d, evaluated at (X,Y), on the 
%    time set t.
%
% inputs:  X,Y are time parameters contained in a [p,q] matrix;
%          repr. is a structure containing 
%            - repr.d has [ndp,m,n] entries (where ndp was previously the
%              number of data points to fit);
%            - repr.a has [3,m,n] entries.
%          t is the time set (this is a [ndp,2] matrix)
%
% outputs: tps is a [p,q,m,n] matrix.
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

  % parameters and defense
  %[mt,nt] = size(X);
  
  a = repr.a;
  d = repr.d;
  
  if isvector(X)
    sX = length(X);
  else
    sX = size(X);
  end
  X = X(:);
  Y = Y(:);
  assert(length(X) == length(Y),'surfaceFitting:dimCheck','X and Y must be of same size');
  
  [ndp,m,n] = size(d);
  assert(3 == size(a,1),'TPSBuild:Asize','a must contain only 3 rows');
  assert(m == size(a,2) && n == size(a,3),'TPSBuild:ADdimCheck','a and d must contain objects of same dimension');
  
  % construction of the solution
  %tps = zeros(mt*nt,m,n);
  tps = zeros(prod(sX),m,n);
  for i = 1:length(X)
    temp = zeros(1,m,n);
    ti = [X(i),Y(i)];
    % radial basis part
    for j = 1:ndp
      r = norm(t(j,:) - ti);
      temp = temp + d(j,:,:).*rbf(r);
    end
    % affine complement
    tps(i,:,:) = temp + a(1,:,:) + a(2,:,:).*X(i) + a(3,:,:)*Y(i);
  end
  tps = reshape(tps,[sX,m,n]);
end

function tangentValues = pointsLifting(M,rootPoints,points,coarse,pointsCoords)
% Lifts the points in the tangent space of the last element of rootPoints 
% by transporting the points through all elements of rootPoints.
%
% function tangentValues = pointsLifting(M,rootPoints,points)
%    returns tangentValues, the representation in T_{rootPoints{end}}M 
%    of the points, lifted thanks to the M.log operator of M, a manifold
%    manopt-structure, and transported through M.transp.
%    The points are first lifted in the first element of rootPoints,
%    and then the tangent values are transported from element to element
%    of rootPoints, and the distance between each rootPoint is
%    accumulated at each transport.
%
% inputs:  M, a manopt manifold-structure containing the M.log operator.
%          rootPoints, a m-cell of point on the manifold M.
%          points, [DIM]-cell of points to be lifted to the tangent 
%            space at rootPoint.
%          coarse (optional) is a boolean that specifies if the points
%            must be approximated by a linear surface within a patch
%            (default: 0).
%            coarse and pointsCoords must be provided together.
%          pointsCoords (optional) are the coordinates of the points to
%            be lifted with the polynomial representation (default: empty).
%            coarseDegree and pointsCoords must be provided together.
%
% outputs: tangentValues, a [DIM]-cell containing the lifted version
%            the points to the tangent spaces at root.
%
% Original author: 
% 	Pierre-Yves Gousenbourger, Jan. 22, 2020.
% Contributors: 
%
% Change log:
% 	Jan. 22, 2020 (PYG) - First version.
%   Aug. 25, 2020 (PYG) - Adaptation to coarseLifting.

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

  % defense
  if nargin <= 3
    coarse = false;
  end
  if ~iscell(rootPoints)
    temp = rootPoints;
    rootPoints = cell(1,1);
    rootPoints{1} = temp;
    clear temp;
  end
  if ~iscell(points)
    temp = points;
    points = cell(1,1);
    points{1} = temp;
    clear temp;
  end
  
  %rootPoints
  assert(isvector(rootPoints),'surfaceFitting:vectorCell','rootPoints must be a scalar or vectorial cell');
  
  % parameter
  sPoints     = size(points);
  sRootPoints = length(rootPoints);
  if sRootPoints > 1
    if isfield(M,'isotransp')
      tp = M.isotransp;
    else
      tp = M.transp;
    end
  end
  
  % method
  points  = points(:);
  nPoints = length(points);
  tangentPoints = cell(size(points));
  coarseGraining = coarse && (nPoints > 2);
  
  % step 1: lift all the dataPoints to obtain tangentPoints in the closest rootPoint.
  for j = 1:nPoints
    tangentPoints{j} = M.log(rootPoints{1},points{j});
  end
  
  % step 2: if necessary, compute the linear surface that best fits the points
  % (only if more than 3 points, of course, otherwise there is no added value)
  if coarseGraining
    tangentPoints = linSurfFit(tangentPoints,pointsCoords);
  end
  
  % step 3: move the coefficients of the polynomial/points to end rootPoints
  tangentValues = cell(size(tangentPoints));
  for j = 1:length(tangentPoints)
    temp = tangentPoints{j};
    for k = 2:length(rootPoints)
      temp = tp(rootPoints{k-1},rootPoints{k},temp) + M.log(rootPoints{k},rootPoints{k-1});
    end
    tangentValues{j} = temp;
    clear temp;
  end
  
  % step 4: if necessary, evaluate the linear surface at the coordinates
  % (only if more than 3 points, of course, again)
  if coarseGraining
    approxValues = cell(size(points));
    for j = 1:length(points)
      approxValues{j} = linSurfVal(tangentValues,pointsCoords(j,:));
    end
    tangentValues = approxValues;
  end
  
  % final reshape
  tangentValues = reshape(tangentValues,sPoints);
end


% Fit a linear surface to the dataPoints
function coefficients = linSurfFit(points,coords)
  sPoint = size(points{1});
  lPoint = prod(sPoint);
  
  % Solver here
  f  = @(x) [lsqe(x,points,coords)];
  % options
  options = optimoptions('fminunc');
  options.OptimalityTolerance = 1e-9;
  options.StepTolerance = 1e-8;
  options.Display = 'off';
  % start point
  x0 = zeros(3*lPoint,1);
  % search
  x  = fminunc(f,x0,options);
  
  % put coefficients in cells
  coefficients{1} = reshape(x(         1:  lPoint),sPoint);
  coefficients{2} = reshape(x(  lPoint+1:2*lPoint),sPoint);
  coefficients{3} = reshape(x(2*lPoint+1:3*lPoint),sPoint);
end

% handle function to compute the least squares error of the fitting surface
function [e,ge] = lsqe(x,points,coords)
  sPoint = size(points{1});
  lPoint = prod(sPoint);
  nData  = length(points);
  
  % transform x into a cell for val use
  range = [         1:  lPoint;
             lPoint+1:2*lPoint;
           2*lPoint+1:3*lPoint];
  coeff{1} = reshape(x(range(1,:)),sPoint);
  coeff{2} = reshape(x(range(2,:)),sPoint);
  coeff{3} = reshape(x(range(3,:)),sPoint);
  
  % initialization
  e = 0;
  ge = zeros(3*lPoint,1);
  
  % evaluation
  for i = 1:nData
    data = points{i};
    val  = linSurfVal(coeff,coords(i,:));
    
    % function e 
    dataNorm = norm(data - val,'fro').^2;
    e = e + dataNorm;
    
    % grad ge
    ge(range(1,:)) = ge(range(1,:)) - 2*(data(:) - val(:));
    ge(range(2,:)) = ge(range(2,:)) - 2*(data(:) - val(:)).*coords(i,1);
    ge(range(3,:)) = ge(range(3,:)) - 2*(data(:) - val(:)).*coords(i,2);
  end
end

% Evaluate the linear surface at the points coordinates.
function val = linSurfVal(coeff,coord)
  val = coeff{1} + coord(1)*coeff{2} + coord(2)*coeff{3};
end

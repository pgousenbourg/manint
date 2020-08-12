function y = karcherMean(M,points,weights)
% returns the 2D-weighted mean of points associated with weights.
% The mean is computed with a so-called Karcher mean.
%
% function Y = karcherMean(M,points)
%    returns Y, the mean of four points with equal weights.
%
% function Y = karcherMean(M,points,weights)
%    returns Y, the mean of four points associated with given weights.
%
% inputs:  M is a manopt manifold-structure containing the exp and log
%            operators.
%          points is a 4-cell of [p,DIM] manifold-valued points in
%            M. Each p sequence of points is associated with weights.
%          weights is a [p,4]-matrix with 4 columns of weights 
%            associated to each point.
%
% outputs: y is the karcher mean of the points wrt the weights, stored 
%          in a [p,DIM]-matrix.
%
% Original author: 
% 	Pierre-Yves Gousenbourger, Dec. 18, 2019.
% Contributors: 
%
% Change log:
%   Jan. 07, 2020 (PYG) - First version based on tensorMean.m

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
  
  % defense on points
  assert(iscell(points),'surfaceFitting:dataNotInCell','Points should be in a cell.');
  assert(length(points) == 4,'blendSurf:fourEntries','There may only be 4 points by patch.');
  cellfun(@(pts) assert(sum(sum(size(points{1}) - size(pts))) == 0,'surfaceFitting:dimCheck','points must all have the same size'),points)
  
  % complete with missing weights
  if nargin == 2
    weights = ones(size(points{1},1),4)./4;
  end
  
  % defense on weights
  assert(ismatrix(weights),'surfaceFitting:dataNotInMatrix','Weights should be in a px4 matrix.');
  assert(size(weights,2) == 4,'blendSurf:fourEntries','There must be 4 columns of weights');
    
  % correct points if only one vector of points
  sW = size(weights,1);
  if sW == 1
    newPoints = cell(1,4);
    for i = 1:4
      temp(1,:,:) = points{i};
      newPoints{i} = temp;
    end
    points = newPoints;
    clear temp newPoints;
  end  
  
  % parameters
  sP = size(points{1});
  p = sP(1);
  dimM = sP(2:end);
  %dimM = M.size();
  if isscalar(dimM)
    dimM = [dimM 1];
  end
  
  % final defense
  assert(size(weights,1) == p,'surfaceFitting:dimCheck','Dimensions mismatch between weights and points');
  
  % Karcher mean with manopt
  y = zeros(sP);
  problem.M = M;
  for i = 1:p
    b = getPoints(points,i,dimM);
    w = getWeights(weights,i);
    
    problem.cost = @(x) karcher(M,x,b,w);
    problem.grad = @(x) dkarcher(M,x,b,w);
    
    options.verbosity = 0;
    warning('off', 'manopt:getHessian:approx') 
    
    y(i,:,:) = trustregions(problem,[],options);
  end
  
  % reshape if necessary
  if sW == 1
    y = reshape(y,dimM);
  end
end


function b = getPoints(points,i,dimM)
  b = cell(1,4);
  for k = 1:4
    temp = points{k};
    b{k} = reshape(squeeze(temp(i,:,:)),dimM);
  end
end

function w = getWeights(weights,i);
  w = weights(i,:);
end

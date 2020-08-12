function y = tensorMean(M,points,weights,varargin)
% returns the tensorized 2D-mean of points associated with weights.
% The mean is computed with diadic averaging.
%
% function Y = tensorMean(M,points)
%    returns Y, the mean of four points with equal weights.
%
% function Y = tensorMean(M,points,weights)
%    returns Y, the mean of four points associated with given weights.
%
% function Y = tensorMean(M,points,weights,'orientation',dim)
%    dim can be 1 or 2, and forces the diadic mean to be done
%    in the dim direction. Default: 2.
%
% inputs:  M is a manopt manifold-structure containing the exp and log
%            operators.
%          points is a 4-cell of [p,DIM] manifold-valued points in
%            M. Each p sequence of points is associated with weights.
%          weights is a [p,4]-matrix with 4 columns of weights 
%            associated to each point.
%
% outputs: y is diadic mean of the points wrt the weights, stored in
%          a [p,DIM]-matrix.
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

  % create the parser
  ip = inputParser();
  addOptional(ip,'orientation',2);

  % parse the parser as a structure.
  parse(ip, varargin{:});
  vars = ip.Results;
  
  % defense on points
  assert(iscell(points),'surfaceFitting:dataNotInCell','Points should be in a cell.');
  assert(length(points) == 4,'blendSurf:fourEntries','There may only be 4 points by patch.');
  cellfun(@(pts) assert(sum(sum(size(points{1}) - size(pts))) == 0,'surfaceFitting:dimCheck','points must all have the same size'),points)
  
  % complete with missing weights
  if nargin - 2*length(varargin) == 2
    weights = ones(size(points{1},1),4)./4;
  end
  
  % defense on weights
  assert(ismatrix(weights),'surfaceFitting:dataNotInMatrix','Weights should be in a px4 matrix.');
  assert(size(weights,2) == 4,'blendSurf:fourEntries','There must be 4 columns of weights');
  assert(vars.orientation == 2 || vars.orientation == 1,'tensorMean:orientation','Orientation must be 1 or 2');
    
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
  
  % preprocessing
  if vars.orientation == 1
    idx = [1 3 2 4];
    points  = points(idx);
    weights = weights(:,idx);
  end
  
  % Dyadic operation
  y = zeros(sP);
  w = dyadicWeights(weights);
  
  tensorOne = dyadAv(M,points{1},points{2},weights(:,1),weights(:,2));
  tensorTwo = dyadAv(M,points{3},points{4},weights(:,3),weights(:,4));
  y         = dyadAv(M,tensorOne,tensorTwo,w(:,1),w(:,2));
  
  if sW == 1
    y = reshape(y,dimM);
  end
end


% dyadic average
function y = dyadAv(M,x1,x2,w1,w2)
  p = length(w1);
  y = zeros(size(x1));
  dim(1) = size(x1,2); dim(2) = size(x1,3);
  for i = 1:p
    if (w1(i) == 0 && w2(i) == 0)
      w = 0;
    else
      w = w2(i)./(w1(i) + w2(i));
    end
    xx = reshape(x1(i,:,:),dim);
    yy = reshape(x2(i,:,:),dim);
    y(i,:,:) = M.exp(xx,w.*M.log(xx,yy));
  end
end

% weights must be given in a 4-cell form.
function wNew = dyadicWeights(w);
  wNew(:,1) = (w(:,1) + w(:,2))./sum(w,2);
  wNew(:,2) = (w(:,3) + w(:,4))./sum(w,2);
end

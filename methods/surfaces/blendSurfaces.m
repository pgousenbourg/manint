function y = blendSurfaces(M,rootPoints,tangentValues,weights,average)
% Blends four tangentValues together. 
%
% function Y = blendSurfaces(M,rootPoints,tangentValues,weights)
%    returns Y, the blended version of the four given tangentValues (in
%    the tangent space T_rootPoints(M)) and with respect to given
%    weights.
%
% function Y = blendSurfaces(M,rootPoints,tangentValues,weights,average)
%    returns Y, the blended version of the tangentValues, with a chosen
%    averaging method, furnished as a handle function.
%    By default, the averaging is done by average = @tensorMean.
%
% inputs:  M, a manopt manifold-structure.
%          rootPoints, a 4-cell containing the rootPoints at which
%            the tangentValues are evaluated.
%          tangentValues, a 4-cell containing the four tangent
%            surfaces in the tangent spaces of the rootPoints.
%            The tangent surfaces in each entries of the cell are stored
%            in a [DIMxmxn] matrix, where [m,n] is the dimention of M.
%          weights, a 4-cell matrix containing the weights used to
%            blend the tangentValues once mapped back to the manifold.
%            The entries of the cell are [DIM]-matrices.
%            If weights is a 4-vector, then the same 4 
%            weights will be used for every tangentValues.
%          average (optional), a handleFunction that averages four
%            points with given weights. Default: tensorMean.
%            The function must be defined as
%            
%            f = @(M,points,weights) average(M,points,weights).
%
% outputs: Y, a matrix [DIMxmxn] with all the blended values in M.
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

  if nargin == 4
    average = @tensorMean;
  end
  
  % parameters
  nRoots = length(tangentValues);
  
  % defense
  assert(iscell(rootPoints),'surfaceFitting:dataNotInCell','rootPoints must be in a cell');
  assert(iscell(tangentValues),'surfaceFitting:dataNotInCell','tangentValues must be in a cell');
  assert(nRoots == 4,'blendSurfaces:fourPoints','There must be four tangentValues to blend.');
  assert(length(rootPoints) == nRoots,'blendSurfaces:fourPoints',...
    'there must be as many rootPoints, tangentValues and weights to blend.');
  if iscell(weights)
    assert(length(weights) == nRoots,'blendSurfaces:fourPoints','there must be four weights for each point to blend');
    for k = 1:nRoots
      if isvector(weights{k})
        cellfun(@(tvi) assert(size(weights{k},1) == size(tvi,1),'surfaceFitting:dimCheck','weights and tangent values must have the same number of points'),tangentValues);
      else
        cellfun(@(tvi) assert(size(weights{k},2) == size(tvi,2) && size(weights{k},1) == size(tvi,1),'surfaceFitting:dimCheck','weights and tangent values must have the same number of points'),tangentValues);
      end
    end
  elseif isvector(weights)
    assert(length(weights) == nRoots,'blendSurfaces:fourPoints','there must be four weights minimum');
  else
    error('invalid weights entry');
  end
  
  % Parameters
  dimTV = size(tangentValues{1});
  [dimGrid,dimM] = getDims(tangentValues{1},rootPoints{1});
  assert(length(dimM) <= 2,'surfaceFitting:tensorManifold','Manifolds must be matrix or vector');
  
  % Transform weights in cell to ensure there is no problem later
  if ~iscell(weights);
    w = cell(size(rootPoints));
    for i = 1:length(w);
      w{i} = repmat(weights(i),dimGrid);
    end
    weights = w;
    clear w;
  end;
  
  % All the actual method is done... here ;-)
  % Mapping back to the manifold
  points = cell(size(rootPoints));
  for k = 1:4
    root = rootPoints{k};
    tVals = reshape(tangentValues{k},[prod(dimGrid),dimM]);
    mVals = zeros(size(tVals));
    for i = 1:size(tVals,1)
      tVali = reshape(tVals(i,:,:),getdimM(dimM));
      mVals(i,:,:) = M.exp(root,tVali);
    end
    points{k} = mVals;
    weights{k} = reshape(weights{k},[prod(dimGrid),1]);
  end
  weights = cell2mat(weights);
  
  % Averaging
  y = average(M,points,weights);
  y = reshape(y,dimTV);
end

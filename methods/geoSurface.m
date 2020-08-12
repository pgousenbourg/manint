function Z = geoSurface(varargin);
% Returns the value Z at coordinate X,Y of the geodesic surface.
%    The surface interpolates a set of data points DATAPOINTS at coordinates 
%    DATACOORDS.
%
%    Z = geoSurface(M,DATAPOINTS,DATACOORDS,X,Y) returns Z, the 
%    value at coordinate [X,Y] of the geodesic surface interpolating 
%    the DATAPOINTS on a manifold M, associated with coordinates DATACOORDS.
%    DATAPOINTS are stored in a vector-cell of N entries, and DATACOORDS
%    are stored in a [Nx2] vector, where each line corresponds to 
%    the coordinate of the associated datapoint.
%
%
% Mandatory entries are:
%             M : The manifold on which the data points are 
%                 expressed.
%
%           X,Y : The grid of coordinates at which the surface
%                 is evaluated.
%
%    dataPoints : The data points to be fitted. Those points
%                 belong to M and are stored in a vector-cell
%                 of N entries.
%
%    dataCoords : Matrix [N x 2] of (x,y)-values, associated
%                 to each dataPoint.
%
%
% Original author: 
% 	Pierre-Yves Gousenbourger, Dec. 19, 2019.
% Contributors: 
%
% Change log:
% 	Dec. 19, 2019 (PYG) - First version.
% 	Mar. 24, 2020 (PYG) - Logic and help.
%   Mar. 25, 2020 (PYG) - Modification for localSurface fitting.
%   Apr.  1, 2020 (PYG) - Modification for geoSurface interpolation.

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

  % ====================================================================
  % input management
  % ====================================================================
  
  % create the parser
  ip = inputParser();
  % required inputs
  addRequired(ip,'M')
  addRequired(ip,'dataPoints')
  addRequired(ip,'dataCoords');
  addRequired(ip,'X')
  addRequired(ip,'Y')
  
  % parse the parser as a structure.
  parse(ip, varargin{:});
  vars = ip.Results;
  
  % Checks inputs, fills empty entries, verifies the coherence of inputs
  % --- dimension tests on dataCoords ----------------------------------
  assert(length(size(vars.dataCoords)) == 2, 'dataCoords must be a matrix');
  assert(isnumeric(vars.dataCoords),'dataCoords must be a matrix of numeric values');
  assert(size(vars.dataCoords,2) == 2, 'dataCoords must have two columns');
  
  % --- dimension tests on dataPoints ----------------------------------
  assert(iscell(vars.dataPoints),'dataPoints must be stored in a cell');
  assert(isvector(vars.dataPoints),'dataPoints must be stored in a vector cell');
  K = size(vars.dataPoints{1});
  for i = 1:length(vars.dataPoints)
  for j = length(K)
    assert(size(vars.dataPoints{i},j) == K(j),'dataPoints must all have the same dimension');
  end
  end
  
  % --- dimension tests on X and Y -------------------------------------
  assert(length(size(vars.X)) == length(size(vars.Y)),'X and Y must have the same number of dimensions');
  for i = 1:length(size(vars.X))
    assert(size(vars.X,i) == size(vars.Y,i),'X and Y must have the same dimensions');
  end
  
  % --- coherence between dataCoords & dataPoints ----------------------
  assert(length(vars.dataPoints) == size(vars.dataCoords,1),'There must be as many dataPoints as dataCoords');
  
  
  % ====================================================================
  % Prepare the parameters
  % ====================================================================
  M          = vars.M;
  dataPoints = vars.dataPoints;
  dataCoords = vars.dataCoords;
  dimX       = size(vars.X);
  dimM       = size(dataPoints{1});
  X          = vars.X(:);
  Y          = vars.Y(:);
  
  % useful variables
  nTimes   = length(X);
  nData    = length(dataPoints);
  domains  = makeDomain(dataCoords);
  
  
  % ====================================================================
  % The algorithm can finally start here
  % ====================================================================
  
  % Generate the geodesics
  Z = zeros([nTimes,dimM]);
  for i = 1:nTimes
    % in which domain are we ?
    domain  = domains(findDomain(X(i),Y(i),domains,dataCoords),2:5);
    data    = dataPoints(domain);
    coords  = dataCoords(domain,:);
    
    % weights
    w = customWeights(X(i),Y(i),coords,@(t) (t));
    
    % geodesics X
    S1 = M.exp(data{1},w(1).*M.log(data{1},data{2}));
    S2 = M.exp(data{3},w(1).*M.log(data{3},data{4}));
     
    % geodesics Y
    Z(i,:,:) = M.exp(S1,w(2).*M.log(S1,S2));
  end
  
  % return values
  Z = squeeze(reshape(Z,[dimX dimM]));
  
end

% === Auxiliary functions ==============================================

% weighting functions :
% g is the 1D function used in the weighting technique
function w = customWeights(xx,yy,rC,g)
  
  % computation of the weights
  t = abs([xx,yy] - rC(1,:));
    
  % width of the windows in X and Y
  width = [rC(2,1) - rC(1,1), rC(3,2) - rC(1,2)];
    
  % defense
  assert(sum(sum((width == 0)))==0,'surfaceWeights:gridMisAlignment','The given rootCoords implied a 0-sized width.\nIs the grid turned by 90 degrees?');
  
  w = g(t./width);
end

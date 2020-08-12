function Z = TSlinearSurface(varargin);
% Returns the value Z at coordinate X,Y of the piecewise linear surface 
% computed in the tangent space of ROOTCOORD.
%    The surface interpolates a set of data points DATAPOINTS at coordinates 
%    DATACOORDS. The tangent space in which the tangent surface is 
%    approximated before being mapped back to M is given by the root 
%    ROOTPOINT.
%
%    Z = TSlinearSurface(M,DATAPOINTS,DATACOORDS,X,Y) returns Z, the 
%    value at coordinate [X,Y] of the piecewise linear surface interpolating 
%    the DATAPOINTS on a manifold M, associated with coordinates DATACOORDS.
%    DATAPOINTS are stored in a vector-cell of N entries, and DATACOORDS
%    are stored in a [Nx2] vector, where each line corresponds to 
%    the coordinate of the associated datapoint. By default, the ROOTPOINT
%    is DATAPOINT{1}.
%
%    Z = TSlinearSurface(M,DATAPOINTS,DATACOORDS,X,Y,'rootPoint',ROOTPOINT)
%    interpolates the datapoints with rootpoint specified at ROOTPOINT.
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
% Optional entries are:
%  rootPoint : the rootPoint in the tangent space of which the local 
%              surface must be computed. The rootPoint is a point on M.
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
%   Apr.  1, 2020 (PYG) - Modification for TSlinearSurface interpolation.

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
  
  % optional inputs
  addOptional(ip,'rootPoint',[]); % default will be checked afterwards

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
  for j = 1:length(K)
    assert(size(vars.dataPoints{i},j) == K(j),'dataPoints must all have the same dimension');
  end
  end
  
  % --- dimension tests on X and Y -------------------------------------
  assert(length(size(vars.X)) == length(size(vars.Y)),'X and Y must have the same number of dimensions');
  for i = 1:length(size(vars.X))
    assert(size(vars.X,i) == size(vars.Y,i),'X and Y must have the same dimensions');
  end
  
  % fill inputs
  if isempty(vars.rootPoint)
    vars.rootPoint = vars.dataPoints{1};
  end
  
  % --- coherence between dataCoords & dataPoints ----------------------
  assert(length(vars.dataPoints) == size(vars.dataCoords,1),'There must be as many dataPoints as dataCoords');
    
  % --- coherence between dataPoints & rootPoints ----------------------
  K = size(vars.dataPoints{1});
  for i = 1:length(K)
    assert(size(vars.rootPoint,i) == K(i),'dataPoints and rootPoint must belong to the same space');
  end
  clear K;
  
  
  
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
  rootPoint  = vars.rootPoint;
  
  % useful variables
  nTimes   = length(X);
  nData    = length(dataPoints);
  
  
  % ====================================================================
  % Generate the methods
  % ====================================================================
  lifting = @(man,root,dp) basicLifting(man,root,dp);
  
  
  % ====================================================================
  % The algorithm can finally start here
  % ====================================================================
  
  % lift the data points to corners of the domain
  tangentPoints = lifting(M,rootPoint,dataPoints);
  
  % generate the fitting curve on T_root M
  tangentSurface = piecewiseLinear(tangentPoints,dataCoords,X,Y);
  
  % Map the solution back to M
  Z = zeros([nTimes,dimM]);
  for i = 1:nTimes
    Z(i,:,:) = M.exp(rootPoint,reshape(tangentSurface(i,:,:),dimM));
  end
  
  % return values
  Z = squeeze(reshape(Z,[dimX dimM]));
end

% === Auxiliary functions ==============================================

% Compute the linear piecewise function
% Output as a [p,q,m,n]-matrix, where (m,n) = dimM, and (p,q) = size(X)
function S = piecewiseLinear(dP,dC,X,Y)
  
  % variables
  if isvector(X)
    sX = length(X);
  else
    sX = size(X);
  end
  [m,n] = size(dP{1});
  
  % coordinates and domains
  X = X(:); Y = Y(:);
  domains = makeDomain(dC);
  
  % method
  S = zeros(prod(sX),m,n);
  for i = 1:prod(sX)
    domain  = domains(findDomain(X(i),Y(i),domains,dC),2:5);
    corners = dP(domain);
    coords  = dC(domain,:);
    weights = customWeights(X(i),Y(i),coords,@(t) (1-t));
    S(i,:,:) = linPlane(corners,weights);
  end
  S = reshape(S,[sX,m,n]);
end


% Compute the linear mean of four points, given the weights w
function y = linPlane(data,w)
  y = zeros(size(data{1}));
  for i = 1:length(data)
    y = y + w(i).*data{i};
  end
end


% weighting functions :
% g is the 1D function used in the weighting technique
function w = customWeights(xx,yy,rC,g)
  
  % parameters  
  sX = size(xx);
  lX = prod(sX);
  
  % computation of the weights
  for i = 1:2;
  for j = 1:2;
    k  = rootPos(i,j);
    t = abs([xx(:),yy(:)] - rC(k,:));
    
    % width of the windows in X and Y
    width = [0 0];
    width(1,1) = rC(rootPos(2,j),1) - rC(rootPos(1,j),1);
    width(1,2) = rC(rootPos(i,2),2) - rC(rootPos(i,1),2);
    width = repmat(width,[lX,1]);
    
    % defense
    assert(sum(sum((width == 0)))==0,'surfaceWeights:gridMisAlignment','The given rootCoords implied a 0-sized width.\nIs the grid turned by 90 degrees?');
    
    % weights for each entry reshaped at the size of X,Y 
    w{k} = reshape(prod(g(t./width),2),sX);
  end
  end
  
  % case where there is only one point
  if length(sX) == 2 && prod(sX) == 1
    w = cell2mat(w);
  end

end

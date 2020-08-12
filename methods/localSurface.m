function [Z, varargout] = localSurface(varargin);
% Returns the value Z of the local surface at coordinate X,Y.
%    The local surface fits a set of data points DATAPOINTS
%    at coordinates DATACOORDS with a fitting focus LAMBDA. The 
%    tangent space in which the tangent surface is approximated before
%    being mapped back to M is given by the root ROOTPOINT.
%
%    Z = localSurface(M,DATAPOINTS,DATACOORDS,X,Y) returns Z, the 
%    value at coordinate [X,Y] of the local surface interpolating the 
%    DATAPOINTS on a manifold M, associated with coordinated DATACOORDS.
%    DATAPOINTS are stored in a vector-cell of N entries, and DATACOORDS
%    are stored in a [Nx2] vector, where each line corresponds to 
%    the coordinate of the associated datapoint. By default, the ROOTPOINT
%    is DATAPOINT{1}.
%
%    Z = localSurface(M,DATAPOINTS,DATACOORDS,X,Y,'rootPoint',ROOTPOINT)
%    interpolates the datapoints with rootpoint specified at ROOTPOINT.
%
%    Z = localSurface(M,DATAPOINTS,DATACOORDS,X,Y,'lambda',LAMBDA)
%    fits the datapoints with a specific focus lambda (from 0 to NaN).
%    When lambda is NaN, interpolation is performed. If lambda is low,
%    then less importance is given to data points.
%
%    Z = localSurface(M,DATAPOINTS,DATACOORDS,X,Y,'options',OPTIONS)
%    permits to add some options to the system. The options are 
%    given hereunder.
%
%    [Z,elements] = localSurface(M,DATAPOINTS,DATACOORDS,X,Y), returns
%    also all the elements used in the blending procedure (dataPoints, 
%    dataCoords, options, lambda, etc.) in a structure ELEMENTS.
%
%    Z = localSurface(M,DATAPOINTS,[],X,Y) tries to interpolate the 
%    datapoints DATAPOINTS without the knowledge of the data coordinates.
%    In that case, the algorithm tries to associate data coordinates to 
%    the data points, on a regular grid. If not possible, an error is 
%    returned.
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
%     lambda : the fitting parameter (scalar), from 0 to NaN, strictly
%              greater than zero. If not provided or NaN, then
%              interpolation is performed by default. 
%
%    options : a structure of option.
%
%
% The possible OPTIONS are the following:
%
%        generationMethod : method that generates the minimal 
%                           representation of the fitting surface on the 
%                           tangent space. Given as a STRING value.
%                           Possible values are: 
%                              'TPS'    - thin plate splines (default)
%                              'Bezier' - bezier surfaces
%
%    reconstructionMethod : method that reconstruct the fitting surface
%                           based on its minimal representation.
%                           Possible values are:
%                              'TPS'    - thin plate splines (default)
%                              'Bezier' - bezier surfaces
%                           /!\ TPS generation does not work with Bezier
%                               reconstruction, and vice-versa.
%
%                 display : boolean value. 1 for displaying stuff, 
%                           0 for shutting it down.
%                           By default: 0 (no display)
%
%                 verbose : Level of verbosity (0 for no verbose, 1 for
%                           normal verbosity or 2 for full verbosity).
%                           By default: 1 (normal verbose)
%                           /!\ Not Yet In Production /!\
%
% Original author: 
% 	Pierre-Yves Gousenbourger, Dec. 19, 2019.
% Contributors: 
%
% Change log:
% 	Dec. 19, 2019 (PYG) - First version.
% 	Mar. 24, 2020 (PYG) - Logic and help.
%   Mar. 25, 2020 (PYG) - Modification for localSurface fitting.

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
  addOptional(ip,'lambda',NaN); % interpolation by default
  addOptional(ip,'options',[]); % default options

  % parse the parser as a structure.
  parse(ip, varargin{:});
  vars = ip.Results;
  
  % Checks inputs, fills empty entries, verifies the coherence of inputs
  checkInputsLocal(vars);
  vars = fillInputsLocal(vars);
  coherenceOfInputsLocal(vars);
  
  
  
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
  lambda     = vars.lambda;
  rootPoint  = vars.rootPoint;
  
  % useful variables
  nTimes   = length(X);
  nData    = length(dataPoints);
  
  
  % ====================================================================
  % Generate the methods
  % ====================================================================
  switch vars.options.generationMethod
    case 'TPS';    generation = @(tp,dc) thinPlateSplineGeneration(tp,dc,lambda);
    case 'bezier'; error('Bezier spline generation is not integrated yet...'); %generation = @(tp,dc) bezierSplineGeneration(tp,dc,lambda);
    otherwise; error('Wrong generation method label provided');
  end
  switch vars.options.reconstructionMethod
    case 'TPS';    reconstruction = @(x,y,ad) thinPlateSplineReconstruction(x,y,ad,dataCoords);
    case 'bezier'; error('Bezier spline reconstruction is not integrated yet...'); %reconstruction = @(x,y,b) bezierSplineReconstruction(x,y,b);
    otherwise; error('Wrong reconstruction method label provided');
  end
  lifting = @(man,root,dp) basicLifting(man,root,dp);
   
  
  % ====================================================================
  % The algorithm can finally start here
  % ====================================================================
  
  % lift the data points to corners of the domain
  tangentPoints = lifting(M,rootPoint,dataPoints);
  
  % generate the fitting curve on T_root M
  TPS = generation(tangentPoints,dataCoords);
  tangentSurface = reconstruction(X,Y,TPS);
  
  % Map the solution back to M
  Z = zeros([nTimes,dimM]);
  for i = 1:nTimes
    Z(i,:,:) = M.exp(rootPoint,reshape(tangentSurface(i,:,:),dimM));
  end
  
  % return values
  Z = squeeze(reshape(Z,[dimX dimM]));
  
  if nargout > 1
    elements = vars;
    elements.options.reconstruction = reconstruction;
    elements.options.generation = generation;
    elements.options = orderfields(elements.options);
    varargout{1} = elements;
  end
end

function vars = fillInputs(vars);
% Fills the inputs for blendedSurfaces.
%
% function vars = fillInputs(vars)
%    Checks all the inputs from the vars and fills empty entries.
%
% The default values if the entries are empty are:
%                       M : returns an error
%                     X,Y : returns an error
%              dataPoints : returns an error
%              dataCoords : returns an error
%                  lambda : NaN (interpolation)
%              rootPoints : takes the dataPoint the closest to 
%                           rootCoords(i);
%              rootCoords : creates a regular grid of the size of 
%                           dataPoints. If not possible, returns an 
%                           error.
%                 domains : generates the domain based on rootCoords, if
%                           possible.
%        generationMethod : 'TPS'
%    reconstructionMethod : 'TPS'
%           liftingMethod : 'basic'
%         blendingWeights : @blendSurfWeights
%          blendingMethod : 'tensor'
%                 display : 0 (no display)
%                 verbose : 1 (normal verbose)
%
% Original author: 
% 	Pierre-Yves Gousenbourger, Mar. 24, 2020.
% Contributors: 
%
% Change log:
% 	Mar. 24, 2020 (PYG) - First version.

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
  

  % --- Fill the lambda ------------------------------------------------
  if isnan(vars.lambda)
    vars.lambda = 0;
  else
    vars.lambda = 1/vars.lambda;
  end
  % --- Fill the rootPoints and rootCoords -----------------------------
  % If none, try to check is dataPoints can suffice.
  % If only rootCoords, use the closest dataPoint as rootPoint.
  % If only rootPoints, return an error.
  if ~isfield(vars.options,'rootPoints') && ~isfield(vars.options,'rootCoords')
    % check if possible to use dataPoints.
    if isRegularGrid(vars.dataCoords)
      vars.options.rootCoords = vars.dataCoords;
      vars.options.rootPoints = vars.dataPoints;
    else
      error('Can''t fallback to dataPoints as rootPoints. Please provide regular grid of rootCoords');
    end
  elseif ~isfield(vars.options,'rootPoints') && isfield(vars.options,'rootCoords')
    % find the closest data point for each rootCoord.
    vars.options.rootPoints = cell(size(vars.options.rootCoords,1),1);
    for i = 1:size(vars.options.rootCoords,1);
      X = vars.options.rootCoords(i,1);
      Y = vars.options.rootCoords(i,2);
      % Find the closest (X,Y).
      distance = sum(([X,Y] - vars.dataCoords).^2,2);
      [~,idx] = min(distance);
      vars.options.rootPoints(i) = vars.dataPoints(idx);
      clear distance;
    end
  elseif ~isfield(vars.options,'rootCoords') && isfield(vars.options,'rootPoints');
    error('Can''t setup coordinates for the provided rootPoints. Please provide rootCoords.');
  end
  
  % --- fill the domains -----------------------------------------------
  if ~isfield(vars.options,'domains');
    vars.options.domains = makeDomain(vars.options.rootCoords);
  end
  % -- reconstruct not given while generation given (or the other way around)
  if isfield(vars.options,'generationMethod') && ~isfield(vars.options,'reconstructionMethod')
    vars.options.reconstructionMethod = vars.options.generationMethod;
  end
  if ~isfield(vars.options,'generationMethod') && isfield(vars.options,'reconstructionMethod')
    vars.options.generationMethod = vars.options.reconstructionMethod;
  end
  % --- fill the generation method on tangent space --------------------
  if ~isfield(vars.options,'generationMethod')
    vars.options.generationMethod = 'TPS';
  end
  % --- fill the reconstruction method on tangent space ----------------
  if ~isfield(vars.options,'reconstructionMethod')
    vars.options.reconstructionMethod = 'TPS';
  end
  % --- fill the lifting method to tangent space -----------------------
  if ~isfield(vars.options,'liftingMethod') && ~isfield(vars.options,'lifting');
    vars.options.liftingMethod = 'basic';
  end
  % --- fill the blending weights computation function -----------------
  if ~isfield(vars.options,'blendingWeights')
    vars.options.blendingWeights = @(tx,ty,rc) blendSurfWeights(tx,ty,rc);
  end
  % --- fill the blending method ---------------------------------------
  if ~isfield(vars.options,'blendingMethod') && ~isfield(vars.options,'blending');
    vars.options.blendingMethod = 'tensor';
  end
  % --- display --------------------------------------------------------
  if ~isfield(vars.options,'display')
    vars.options.display = 0; % no display by default
  end
  % --- verbose --------------------------------------------------------
  if ~isfield(vars.options,'verbose')
    vars.options.verbose = 1; % normal verbose by default
  end
end

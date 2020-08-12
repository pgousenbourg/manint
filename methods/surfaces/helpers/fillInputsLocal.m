function vars = fillInputsLocal(vars);
% Fills the inputs for localSurfaces.
%
% function vars = fillInputsLocal(vars)
%    Checks all the inputs from the vars and fills empty entries.
%
% The default values if the entries are empty are:
%                       M : returns an error
%                     X,Y : returns an error
%              dataPoints : returns an error
%              dataCoords : returns an error
%                  lambda : NaN (interpolation)
%               rootPoint : takes the first dataPoint;
%        generationMethod : 'TPS'
%    reconstructionMethod : 'TPS'
%                 display : 0 (no display)
%                 verbose : 1 (normal verbose)
%
% Original author: 
% 	Pierre-Yves Gousenbourger, Mar. 24, 2020.
% Contributors: 
%
% Change log:
% 	Mar. 24, 2020 (PYG) - First version.
%   Mar. 25, 2020 (PYG) - Modification for localSurface

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
  % --- Fill the rootPoint ---------------------------------------------
  if isempty(vars.rootPoint)
    vars.rootPoint = vars.dataPoints{1};
  end
  
  % --- fill the generation method on tangent space --------------------
  if ~isfield(vars.options,'generationMethod')
    vars.options.generationMethod = 'TPS';
  end
  % --- fill the reconstruction method on tangent space ----------------
  if ~isfield(vars.options,'reconstructionMethod')
    vars.options.reconstructionMethod = 'TPS';
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

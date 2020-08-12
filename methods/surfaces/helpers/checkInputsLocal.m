function checkInputsLocal(vars);
% Checks the consistency of the inputs for localSurface.
%
% function checkInputsLocal(vars)
%    Checks all the inputs from the vars, independently.
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
  
  % --- tests on lambda ------------------------------------------------
  assert(vars.lambda > 0 || isnan(vars.lambda),'lambda must be NaN or strictly greater than zero.');
  
  % --- name of the generation method ----------------------------------
  if isfield(vars.options,'generationMethod')
    str = vars.options.generationMethod;
    assert(strcmp(str,'TPS') || strcmp(str,'bezier'),'incorrect value for generationMethod');
  end
  
  % --- name of the reconstruction method ------------------------------
  if isfield(vars.options,'reconstructionMethod')
    str = vars.options.reconstructionMethod;
    assert(strcmp(str,'TPS') || strcmp(str,'bezier'),'incorrect value for reconstructionMethod');
  end
  
  clear str;
  
  % --- display --------------------------------------------------------
  if isfield(vars.options,'display')
    assert(vars.options.display == 0 || vars.options.display == 1,'incorrect value for display');
  end
  
  % --- verbose --------------------------------------------------------
  if isfield(vars.options,'verbose')
    verb = vars.options.verbose;
    assert(verb == 0 || verb == 1 || verb == 2,'incorrect value for verbose');
  end
end

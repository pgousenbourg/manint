function checkInputs(vars);
% Checks the consistency of the inputs for blendedSurfaces.
%
% function checkInputs(vars)
%    Checks all the inputs from the vars, independently.
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
  
  % --- dimension tests on rootPoints ----------------------------------
  if isfield(vars.options,'rootPoints') 
    assert(iscell(vars.options.rootPoints),'rootPoints must be stored in a cell');
    assert(isvector(vars.options.rootPoints),'rootPoints must be stored in a vector cell');
    K = size(vars.options.rootPoints{1});
    for i = 1:length(vars.options.rootPoints)
    for j = length(K)
      assert(size(vars.options.rootPoints{i},j) == K(j),'rootPoints must all have the same dimension');
    end
    end
  end
  
  % --- dimension tests on rootCoords ----------------------------------
  if isfield(vars.options,'rootCoords')
    assert(length(size(vars.options.rootCoords)) == 2, 'rootCoords must be a matrix');
    assert(isnumeric(vars.options.rootCoords),'rootCoords must be a matrix of numeric values');
    assert(size(vars.options.rootCoords,2) == 2, 'rootCoords must have two columns');
    assert(isRegularGrid(vars.options.rootCoords));
  end
  
  % --- dimension tests on domains -------------------------------------
  if isfield(vars.options,'domains');
    assert(length(size(vars.options.domains)) == 2, 'domains must be a matrix');
    assert(isnumeric(vars.options.domains),'domains must be a matrix of numeric values');
    assert(size(vars.options.domains,2) == 5, 'domains must have five columns');
  end
  
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
  
  % --- name of the lifting method -------------------------------------
  if isfield(vars.options,'liftingMethod')
    str = vars.options.liftingMethod;
    assert(strcmp(str,'basic') || strcmp(str,'transport') || strcmp(str,'coarse'),'incorrect value for liftingMethod');
  end
  
  % --- name of the blending method ------------------------------------
  if isfield(vars.options,'blendingMethod')
    str = vars.options.blendingMethod;
    assert(strcmp(str,'tensor') || strcmp(str,'karcher'),'incorrect value for blendingMethod');
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

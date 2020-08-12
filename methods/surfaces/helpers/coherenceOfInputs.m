function coherenceOfInputs(vars);
% Checks the consistency of the inputs for blendedSurfaces.
%
% function coherenceOfInputs(vars)
%    Checks the coherence between pairs of inputs from the vars.
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

    
  % --- coherence between dataCoords & dataPoints ----------------------
  assert(length(vars.dataPoints) == size(vars.dataCoords,1),'There must be as many dataPoints as dataCoords');
  
  % --- coherence between rootCoords & rootPoints ----------------------
  assert(length(vars.options.rootPoints) == size(vars.options.rootCoords,1),'There must be as many rootPoints as rootCoords');
  
  % --- coherence between dataPoints & rootPoints ----------------------
  K = size(vars.dataPoints{1});
  for i = 1:length(K)
    assert(size(vars.options.rootPoints{1},i) == K(i),'dataPoints and rootPoints must belong to the same space');
  end
  
  clear K;
  
  % --- coherence between X,Y and dataCoords, rootCoords ---------------
  minX = min(vars.options.rootCoords(:,1)); maxX = max(vars.options.rootCoords(:,1));
  minY = min(vars.options.rootCoords(:,2)); maxY = max(vars.options.rootCoords(:,2));
  X = vars.X(:); Y = vars.Y(:);
  assert(min(X) >= minX && max(X) <= maxX && min(Y) >= minY && max(Y) <= maxY,'X and Y must be in the domain set by rootCoords');
  
  clear minX minY maxX maxY X Y;
  
  % --- coherence between reconstructionMethod and generationMethod ----
  assert(strcmp(vars.options.reconstructionMethod,vars.options.generationMethod),'generation and reconstruction should be compatible');
  
  % --- coherence between domains and rootCoords -----------------------
  doms = vars.options.domains(2:5);
  doms = doms(:);
  rC   = length(vars.options.rootCoords);
  assert(max(doms) <= rC && min(doms) >= 0, 'domains indexes must correspond to existing rootCoords');
  clear rC doms;
end

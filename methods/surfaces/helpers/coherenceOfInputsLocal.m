function coherenceOfInputsLocal(vars);
% Checks the consistency of the inputs for localSurface.
%
% function coherenceOfInputsLocal(vars)
%    Checks the coherence between pairs of inputs from the vars.
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

    
  % --- coherence between dataCoords & dataPoints ----------------------
  assert(length(vars.dataPoints) == size(vars.dataCoords,1),'There must be as many dataPoints as dataCoords');
    
  % --- coherence between dataPoints & rootPoints ----------------------
  K = size(vars.dataPoints{1});
  for i = 1:length(K)
    assert(size(vars.rootPoint,i) == K(i),'dataPoints and rootPoint must belong to the same space');
  end
  
  clear K;
  
  % --- coherence between reconstructionMethod and generationMethod ----
  assert(strcmp(vars.options.reconstructionMethod,vars.options.generationMethod),'generation and reconstruction should be compatible');
  
end

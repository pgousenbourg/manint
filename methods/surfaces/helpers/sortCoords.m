function [dC,varargout] = sortCoords(dataCoords,column)
% Sorts the dataCoords along the given column in ascending order.
%
%  DC = sortCoords(ROOTCOORDS)
%    returns the list of coordinates with Xs ordered by ascending order.
%
%  DC = sortCoords(ROOTCOORDS,COLUMN)
%    returns the list of coordinates ordered with respect to the given
%    COLUMN (1 or 2). By default, the first column is chosen.
%
%  [DC,I] = sortCoords(ROOTCOORDS)
%    returns also the list of indexes (I) of the sorting.
%
% Inputs:  ROOTCOORDS is a [px2] vector of coordinates;
%          COLUMN is 1 or 2;
%
% Output:  DC is the dataCoords ordered after the 1st or 2nd column.
%
% Original author: 
% 	Pierre-Yves Gousenbourger, Jan. 08, 2020.
% Contributors: 
%
% Change log:
% 	Jan. 08, 2020 (PYG) - First version.

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

  if nargin == 1
    column = 1;
  end
  assert(size(dataCoords,2) == 2,'sortCoords:dimCheck','rootCoords must have two columns');
  assert(column == 1 || column == 2,'sortCoords:column','column must be 1 or 2');
  
  [~,I] = sort(dataCoords(:,column));
  dC = dataCoords(I,:);
  
  if nargout == 2
    varargout{1} = I;
  end
end

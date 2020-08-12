function k = rootPos(i,j)
% Returns the line-position of a matrix-position representation.
%
% function k = rootPos(i,j)
%   returns k, the position in a line-vector when this position is
%   given in (x,y)-coordinates, stored in a p-vector.
%
% function k = rootPos(i)
%   returns k, the matrix-position when the value is given as a
%   line-position i, stored in a [px2]-matrix.
%
% In other words, this helper is doing the bijection:
%         12  ---  22     <->     3  ---  4
%          |       |              |       |
%  y ^     |       |              |       |
%    |    11  ---  21             1  ---  2
%    --> x
%
% Original author: 
% 	Pierre-Yves Gousenbourger, Dec. 18, 2019.
% Contributors: 
%
% Change log:
% 	Dec. 18, 2019 (PYG) - First version.

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

  % defense
  if nargin == 1
    one = (i == 1);
    two = (i == 2);
    three = (i == 3);
    four = (i == 4);
    passed = (sum(one + two + three + four - ones(size(one))) == 0);
    assert(passed,'rootPos:inputType','i must be integer between 1 and 4');
    clear passed one two three four;
  else
    one = (i == 1);
    two = (i == 2);
    passed = (sum(one + two - ones(size(i))) == 0);
    assert(passed,'rootPos:inputType','i must be 1 or 2');
    one = (j == 1);
    two = (j == 2);
    passed = (sum(one + two - ones(size(j))) == 0);
    assert(passed,'rootPos:inputType','j must be 1 or 2');
    clear passed one two;
  end
  
  % parameters
  dimi = size(i);
  i = i(:);
  if nargin == 2
    j = j(:);
  end
  
  % Translation and output
  if nargin == 2
    k = i+(j-1).*2;
  else
    k(:,1) = abs(mod(i,2)-2);
    k(:,2) = ceil(i/2);
  end
  
end

function [cost,route] = dijkstra(graph,in,out)
% Solves the shortest path problem given a graph (stored as an adjacency matrix)
% and in and out nodes (given as indexes).
%
% function [COST,ROUTE] = dijkstra(GRAPH)
%    Returns the min COST and the min-cost ROUTE between all vertices
%    of the graph GRAPH, using Dijkstra algorithm. GRAPH is given as an
%    adjacency matrix.
%
% function [COST,ROUTE] = dijkstra(GRAPH,IN,OUT)
%    Returns the min COST and the min-cost ROUTE between IN and OUT in
%    the given GRAPH, using Dijkstra algorithm. GRAPH is given as an
%    adjacency matrix.
%
% Inputs:  GRAPH, N by N adjacency matrix containing the cost of each
%            edge.
%          IN, OUT (optional), the indexes of the input and output
%            points of the graph. Scalar or vector of length, only.
%            By default, 1:N is given to IN and OUT.
%
% Outputs: COST, a LXM matrix containing the costs of each path for the
%            L inputs and M outputs.
%          ROUTE, a LXM cell containing the paths between the L inputs
%            and the M outputs.
%
% Original author: 
% 	Pierre-Yves Gousenbourger, Jan. 16, 2020.
% Contributors: 
%
% Change log:
% 	Jan. 16, 2020 (PYG) - First version.
%   Jan. 21, 2020 (PYG) - Courtesy to Joseph Kirk (jdkirk630@gmail.com)
%                         who provided the helper function on which this
%                         algorithm is based on (Feb. 02, 2015).

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

  % preprocessing
  n = size(graph,1);
  if nargin == 1
    in = 1:n;
    out = 1:n;
  end
  
  % defense
  assert(size(graph,1) == size(graph,2),'dijkstra:square','The adjacency matrix must be squared.');
  assert(length(in) == prod(size(in)),'dijkstra:dimCheck','Input vector of points must be 1D');
  assert(length(out) == prod(size(out)),'dijkstra:dimCheck','Output vector of points must be 1D');
  assert(max(in) <= n && min(in) >= 1 && max(out) <= n && min(out) >= 1,'dijkstra:inOut','indexes must be between 1 and number of nodes');
  
  % routine
  A = (graph ~= 0);
  [cost,route] = dijkstraMinCost(A,graph,in,out);

end

function graph = domain2adjacency(varargin)
% Transforms a domain to an adjacency matrix of a graph, where the
% cost of each edge is the distance between each connected roots.
%
% function graph = domain2adjacency(M,domains,rootPoints)
%    returns an adjacency matrix between each root of the domain. The
%    cost of each edge is given by the Riemannian distance between
%    the connected rootPoints.
%
% function graph = domain2adjacency(domains,rootCoords)
%    returns an adjacency matrix between each root of the domain. The
%    cost of each edge is given by the Euclidean distance between each
%    node in the domain.
%
% Inputs:  M, a manifold structure from Manopt.
%          DOMAINS, a [px5] matrix of domains (simplified or not) where
%            each line correspond to the four corners of a patch.
%          ROOTPOINTS, a cell with all the roots associated with the domain,
%            expressed on M. There may be n elements.
%          ROOTCOORDS, a matrix [nx2] containing the coordinates in the
%            domain.
%
% Outputs: graph, an adjacency matrix [n x n] with distance costs
%          between each edge.
%
% Original author: 
% 	Pierre-Yves Gousenbourger, Jan. 16, 2020.
% Contributors: 
%
% Change log:
% 	Jan. 16, 2020 (PYG) - First version.
% 	Jan. 22, 2020 (PYG) - Two possible uses of domain2adjacency.

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
  
  % create the parser
  ip = inputParser();
  if nargin == 3
    addRequired(ip,'M')
  end
  addRequired(ip,'domains')
  addRequired(ip,'roots')
  % parse
  parse(ip, varargin{:});
  vars = ip.Results;
  % variables  
  if nargin == 3
    M = vars.M;
    distance = @(x,y) M.dist(x,y);
  elseif nargin == 2
    distance = @(x,y) norm(x-y);
  else
    error('Wrong number of input arguments')
  end
  domains = vars.domains;
  roots = vars.roots;
  
  % defense
  if nargin == 3
    assert(isvector(roots) && iscell(roots),'domains2adjacency:isVector','Roots must be stored in a vectorial cell');
    n = length(roots);
  elseif nargin == 2
    assert(size(roots,2) == 2,'domains2adjacency:isVector','Root coords must be stored in matrix with 2 columns.');
    n = size(roots,1);
    temp = cell(n,1);
    for i = 1:n
      temp{i} = roots(i,:);
    end
    roots = temp; clear temp;
  end
  assert(max(max(domains)) <= n,'domains2adjacency:dimCheck','domains refer to inexistant rootPoints');
  
  % Simplify the domain
  domains = simplifyDomain(domains);
  domains = domains(:,2:5);
  
  % Create the adjacency matrix
  graph = zeros(n,n);
  
  for i = 1:size(domains,1)
    for k = 1:3
    for l = k+1:4
      % indexes of the patch
      in  = domains(i,k);
      out = domains(i,l);
      
      % distance
      dist = distance(roots{in},roots{out});
      
      % storage as a non-oriented graph
      graph(in,out) = dist;
      graph(out,in) = dist;
    end      
    end
  end
end

function plotDomain(domains,rootCoords)
% plots the domain structured by DOMAINS and ROOTCOORDS.
%
% function plotDomain(DOMAINS,ROOTCOORDS)
%    Plots a 2D-plot with the X-Y domains decomposed in subdomains.
%
% See also: makeDomain, refineDomain.
%
% Original author: 
% 	Pierre-Yves Gousenbourger, Jan. 13, 2020.
% Contributors: 
%
% Change log:
% 	Jan. 13, 2020 (PYG) - First version.

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
  assert(size(domains,2) == 5,'findDomain:dimCheck','domains must have five columns');
  assert(size(rootCoords,2) == 2,'findDomain:dimCheck','rootCoords must have two columns');
  assert(length(size(domains)) == 2 && length(size(rootCoords)) == 2,'findDomain:dimCheck','rootCoords and domains must be matrices');
  assert(max(max(domains)) <= size(rootCoords,1) && min(min(domains(:,2:5))) >= 1,'plotDomain:coherence','domains maps inconsistent rootCoords');
  
  hold on;
  title('Plot of the domain');
  for i = 1:size(domains,1);
    limits = computeLimits(i,domains,rootCoords);
    ymin = limits(1);
    xmin = limits(2);
    xmax = limits(3);
    ymax = limits(4);
    toplot = [xmin,ymin;
              xmax,ymin;
              xmax,ymax;
              xmin,ymax;
              xmin,ymin];
    if isSemi(i,domains,rootCoords)
      plot(toplot(:,1),toplot(:,2),':b');
    else
      plot(toplot(:,1),toplot(:,2),'b');
    end
    plot(rootCoords(:,1),rootCoords(:,2),'.b','MarkerSize',20);
  end
  axis equal off;
  hold off
end

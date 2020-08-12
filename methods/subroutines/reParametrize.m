function tOut = reParametrize(t,dataCoords)
% REPARAMETRIZE changes the vector t according to the dataCoords.
%
%    TOUT = REPARAMETRIZE(T,DATACOORDS) reparametrizes T such that the
%    DATACOORDS correspond to integer times i = 0,1,2,... DataCoords and
%    TOUT must be 1D-vectors.
%
%    XOUT = REPARAMETRIZE(X,DATACOORDS) reparametrize X the 2-matrix 
%    X such that DATACOORS coorespond to interger times as well. 
%    DATACOORDS must be a vector of coordinates.
%
%    EXAMPLES:
%      t = [1,2,3,4,5] and dataCoords = [1,3,5];
%      Then, tOut = reParametrize(t,dataCoords) will return
%      tOut = 
%         [0,0.5,1,1.5,2].
%
%      X = [1,2,3,4,5;   and dataCoords = [1,2];
%           1,2,3,4,5]
%      Then Xout = reParametrize(X,dataCoords) will reparametrize as
%      Xout = [0,0.25,0.5,0.75,1;
%              0,0.25,0.5,0.75,1].
%      
% Original author: 
%   Pierre-Yves Gousenbourger, Feb. 27, 2018
% Contributors: 
%
% Change log:
%   Oct. 3, 2019 (PYG):
%      Generalization to 2D.
%

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

    % preprocessing flatten the matrix if needed
    sT = size(t); 
    t = t(:);
    
    % reshaping process
    p = length(t);
    tt = zeros(sT);
    for i = 1:p
      % detect the closest data points
      idx = find(dataCoords - t(i) >= 0);
      if isempty(find(dataCoords == t(i))) % data coord not at t(i)
        t1 = dataCoords(idx(1)-1);
        t2 = dataCoords(idx(1));
        tt(i) = (t(i) - t1)./(t2 - t1) + (idx(1)-2);
      else
        tt(i) = idx(1)-1;
      end
    end
    
    % postprocessing if needed, reshape tt at the matrix size
    tOut = reshape(tt,sT);
end

function [y,varargout] = interpGeo(varargin)
% INTERPGEO returns the value Y of the piecewise geodesic interpolating
% spline at time T.
%    The piecewise geodesic spline interpolates a set of data points 
%    DATAPOINTS at time-parameter DATACOORDS; the spline is composed of
%    geodesics continuously connected at the data points.
%
%    Y = INTERPGEO(M,DATAPOINTS,DATACOORDS,T) returns the
%    value of the piecewise geodesic spline at time T, interpolating the 
%    DATAPOINTS at coordinates DATACOORDS.
%
%    Y = INTERPGEO(M,DATAPOINTS,[],T) returns Y, the value at time T of the 
%    piecewise geodesic spline interpolating the DATAPOINTS on a 
%    manifold M. By default, the method assumes that the DATAPOINTS(i) 
%    is associated with the time parameter DATACOORDS(i) = i.
%
%    [Y,TSPAN] = INTERPGEO(M,DATAPOINTS,[],T) returns also the time needed
%    to compute each value of the fitting curve.
%   
%    Inputs: 
%      M is a manifold structure from Manopt - www.manopt.org.
%      DATAPOINTS must be cell of n entries, where n is 
%        the number of data points.
%        It can also be a matrix [n_rows,n_cols,n], where (n_rows,n_cols) 
%        is the dimension of the search space, and n is the number of data 
%        points.
%      DATACOORDS must be a vector of size d, containing the coordinate 
%        at which each DATAPOINTS must be fitted.
%        Default value: [0 1 2 ... d-1].
%      T can be a scalar or a vector of size p.
%
%    Outputs:
%      Y is a matrix of size [n_rows,n_cols] containing the values of the 
%        interpolating curve evaluated at the p times contained in T.
%      V is the matrix of size [n_rows,n_cols,p] containing the values of
%        the velocity of the curve at T.
%      TSPAN (optional) is a vector of length p that contains the time spent
%        to compute the interpolating curve at times contained in T.
%
% Original author: 
%   Pierre-Yves Gousenbourger, Feb. 27, 2018
% Contributors: 
%
% Change log:
% 	Sep. 28, 2018 (PYG):
%      First version
%   May. 12, 2020 (PYG):
%      DATACOORDS is mandatory

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
addRequired(ip,'M')
addRequired(ip,'dataPoints')
addRequired(ip,'dataCoords')
addRequired(ip,'t')
addOptional(ip,'display',false); % no display by default

% parse the parser as a structure.
parse(ip, varargin{:});
vars = ip.Results;

% defense code
if ~isempty(vars.dataCoords); assert(min(vars.t) >= min(vars.dataCoords) && max(vars.t) <= max(vars.dataCoords),'Please choose t inside the range of dataCoords.'); end
% transform cell to matrix
if iscell(vars.dataPoints)
  assert(isvector(vars.dataPoints),'The data must be in a vector-cell');
  [m,n] = size(vars.dataPoints{1,1});
  dP = zeros(m,n,length(vars.dataPoints));
  for i = 1:length(vars.dataPoints)
    dP(:,:,i) = vars.dataPoints{i};
  end
  vars.dataPoints = dP;
end
clear m n dP
assert(length(size(vars.dataPoints)) == 3, 'dataPoints must be stored in a 3-matrix.');

% variables
M = vars.M;
t = vars.t;
dataPoints = vars.dataPoints;
[n_rows,n_cols,n] = size(dataPoints);
p 		= length(t);
if isempty(vars.dataCoords); 
	dataCoords = linspace(0,n-1,n);
else
	dataCoords = vars.dataCoords;
end

% info display
if vars.display 
	disp('=====================================================');
	fprintf('Piecewise Geodesic spline tool: %26s\n',problem.M.name());
	disp('-----------------------------------------------------');
	disp('Curve computation...');
end

% we should compare the logs beforehand (this can be done offline)
logs = zeros(n_rows,n_cols,n-1);
for i = 1:n-1
    logs(:,:,i) = M.log(dataPoints(:,:,i),dataPoints(:,:,i+1));
end
    

% evaluation at t of the curve
y = zeros(p,n_rows,n_cols);
tSpan = zeros(1,p);

for i = 1:p
  tStart = tic;
	% detect the closest data points
	idx = find(dataCoords - t(i) >= 0);
	if isempty(find(dataCoords == t(i))) % data coord not at t(i)
		a = dataPoints(:,:,idx(1)-1);
		t1 = dataCoords(idx(1)-1);
		t2 = dataCoords(idx(1));
		tt = (t(i) - t1)./(t2 - t1);
		% geodesic 
		y(i,:,:) = M.exp(a,logs(:,:,idx(1)-1),tt);%;
	else
		y(i,:,:) = dataPoints(:,:,idx(1));
  end
  tSpan(i) = toc(tStart);
end

if vars.display;
	fprintf('Piecewise geodesic spline evaluated in: '); 
	fprintf('%17.2f [s] \n\n',sum(tSpan));
end

% Output tSpan if requested.
nout = max(nargout,1) - 1;
if nout == 1
	varargout{1} = curveVelocity(M,y,t);
    varargout{2} = tSpan;
elseif nout == 2
    varargout{2} = tSpan;
end

end

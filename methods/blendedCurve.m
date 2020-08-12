function [y, varargout] = blendedCurve(varargin);
% BLENDEDCURVE returns the value Y of the blended Bezier spline at time T.
%    The blended Bezier spline fits a set of data points DATAPOINTS
%    at time-parameter DATACOORDS with a fitting focus LAMBDA.
%
%    Y = BLENDEDCURVE(M,DATAPOINTS,DATACOORDS,T) returns the
%    value of the blended Bezier spline at time T, interpolating the 
%    DATAPOINTS at coordinates DATACOORDS.
%
%    Y = BLENDEDCURVE(M,DATAPOINTS,[],T). In that case, the method
%    that the DATAPOINTS(i) is associated with the time parameter 
%    DATACOORDS(i) = i-1.
%
%    Y = BLENDEDCURVE(M,DATAPOINTS,DATACOORDS,T,'lambda',LAMBDA) returns 
%    the value of the blended Bezier spline fitting the DATAPOINTS at 
%    given DATACOORDS, but with fitting controlled by the parameter 
%    LAMBDA. It indicates the importance of fitting. A high lambda tends 
%    to interpolation.
%   
%    [Y,TSPAN] = BLENDEDCURVE(M,DATAPOINTS,DATACOORDS,T) returns also 
%    the time needed to compute each value of the fitting curve.
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
%      OPTIONS that can be:
%         'lambda' with a scalar LAMBDA setting the importance to
%            fitting. The highest, the closest to the data points.
%            Default value: 100.
%
%    Outputs:
%      Y is a matrix of size [n_rows,n_cols,n] containing the values of the 
%        blended Bezier spline evaluated at the p times contained in T.
%      V is the matrix of size [n_rows,n_cols,p] containing the values of
%        the velocity of the curve at T.
%      TSPAN (optional) is a vector that contains the time spent to compute
%        the fitting curve at times contained in T.
%
% Original author: 
%   Pierre-Yves Gousenbourger, Feb. 27, 2018
% Contributors: 
%
% Change log:
% 	Sep. 28, 2018 (PYG):
%      First version
%   Nov. 20, 2018 (EM): 
%      Homogenisation of the code w.r.t. the other codes of the repository
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
addOptional(ip,'lambda',[]); % empty by default
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

% add the tools for blending in the path
% addpath(genpath([pwd,'/blended']));

% create the problem structure
pb.M = vars.M;
pb.space = pb.M.name();
pb.data = vars.dataPoints;
t = vars.t;
p = length(t);
[n_rows,n_cols,n] = size(pb.data); % for this manifold, m >> n

% info display
if vars.display; 
	disp('=====================================================');
	fprintf('Bezier blended spline tool: %26s\n',problem.M.name());
	disp('-----------------------------------------------------');
	tStart = tic;
	disp('Curve computation...');
end

% control points generation
if isempty(vars.dataCoords) && isempty(vars.lambda)
	% interpolation (default)
	[pb.control,pb.Tcp] = cp_interp(pb);
elseif isempty(vars.dataCoords) && ~isempty(vars.lambda)
	% fitting
	[pb.control,pb.Tcp] = cp_approx(pb,vars.lambda);
elseif isempty(vars.lambda)
	% interpolation at given dataCoords
	warning('cp_regression not the smartest launch. Launching cp_interp with reparametrization, instead.');
	[pb.control,pb.Tcp] = cp_interp(pb);
	% reparametrization
	t = reParametrize(t,vars.dataCoords);
else
	% fitting at given dataCoords
	warning('cp_regression not the smartest launch. Launching cp_approx with lambda = 10^8 with reparametrization, instead.');
	[pb.control,pb.Tcp] = cp_approx(pb,vars.lambda);
	% reparametrization
	t = reParametrize(t,vars.dataCoords);
end

if vars.display;
	fprintf('Blended spline computed in: ');  
	fprintf('%21.2f [s]\n\n',toc(tStart));
	disp('Curve evaluation...');
end

y = zeros(p,n_rows, n_cols);
tSpan = zeros(1, p);

% evaluation at t of the curve
for i = 1:p
    tStart = tic;
    y(i,:,:) = blended_spline(pb,t(i),vars.display);
    tSpan(i) = toc(tStart);
end

if vars.display;
	fprintf('Blended spline evaluated in: '); 
	fprintf('%17.2f [s] \n\n',tSpan);
end

% Output tSpan if requested.
nout = max(nargout,1) - 1;
if nout == 1
  varargout{1} = curveVelocity(pb.M,y,t);
elseif nout == 2
  varargout{1} = curveVelocity(pb.M,y,t);
  varargout{2} = tSpan;
end

end


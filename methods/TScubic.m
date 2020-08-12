function [y, varargout] = TScubic(varargin);
% TSCUBIC returns the value Y of the cubic spline computed at a given 
% tangent space and mapped back to M (the TS cubic spline) at time T.
%    The spline fits a set of data points DATAPOINTS at time-parameter
%    DATACOORDS with a fitting focus LAMBDA.
%
%    Y = TSCUBIC(M,DATAPOINTS,DATACOORDS,T) returns the value of the 
%    local cubic spline at time T, interpolating the DATAPOINTS at coordinates 
%    DATACOORDS. The ROOTPOINT of the tangent space is by default the 
%    mid-point of the data set.
%
%    Y = TSCUBIC(M,DATAPOINTS,[],T) returns Y, the value at time T of the 
%    TS cubic spline interpolating the DATAPOINTS on a manifold M. By 
%    default, the method assumes that the DATAPOINTS(i) is associated 
%    with the time parameter DATACOORDS(i) = i-1. The ROOTPOINT of the
%    tangent space is by default the mid-point of the data set.
%
%    Y = TSCUBIC(M,DATAPOINTS,[],T,'lambda',LAMBDA) returns the value of
%    the TS cubic spline fitting the DATAPOINTS. LAMBDA indicates 
%    the importance of fitting. A high lambda tends to interpolation.
%    
%    Y = TSCUBIC(M,DATAPOINTS,[],T,'rootPoint',ROOTPOINT) returns the
%    value of the TS cubic spline computed at a specific ROOTPOINT.
%    
%    [Y,TSPAN] = TSCUBIC(M,DATAPOINTS,[],T) returns also the time needed
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
%      OPTIONS that can be:
%         'lambda' with a scalar LAMBDA setting the importance to
%            fitting. The highest, the closest to the data points.
%            Default value: 100.
%         'rootPoint' will set the rootpoint at which the tangent space
%            is considered.
%
%    Outputs:
%      Y is a matrix of size [n_rows,n_cols,p] containing the values of
%        the TS cubic spline evaluated at the p times contained in T.
%      V is the matrix of size [n_rows,n_cols,p] containing the values of
%        the velocity of the curve at T.
%      TSPAN (optional) is a vector that contains the time spent to 
%        compute the fitting curve at times contained in T.
%
% Original author: 
%   Pierre-Yves Gousenbourger, Nov. 21, 2018
% Contributors: 
%
% Change log:
% 	Nov. 21, 2018 (PYG):
%      First version by updating 'blended.m'
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
addOptional(ip,'rootPoint',[]);

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
M = vars.M;
space = M.name();
data = vars.dataPoints;
t = vars.t;
p = length(t);
[n_rows,n_cols,n] = size(data); % for this manifold, m >> n

% Default rootPoint
if ~isempty(vars.rootPoint)
	rootPoint = vars.rootPoint;
else
	rootPoint = data(:,:,floor(n/2));
end

% info display
if vars.display; 
	disp('=====================================================');
	fprintf('TS cubic spline tool: %26s\n',space);
	disp('-----------------------------------------------------');
	tStart = tic;
	disp('Curve computation...');
end

% Creation of the problem structure
pb.M    = euclideanfactory(n_rows,n_cols);
pb.data = zeros(n_rows,n_cols,n);

% Points are lifted to the Tangent Space at rootPoint
for i = 1:n
	pb.data(:,:,i) = M.log(rootPoint,data(:,:,i));
end

% control points generation
if isempty(vars.dataCoords) && isempty(vars.lambda)
	% interpolation (default)
	[controlPoints,~] = cp_interp(pb);
elseif isempty(vars.dataCoords) && ~isempty(vars.lambda)
	% fitting
	[controlPoints,~] = cp_approx(pb,vars.lambda);
elseif isempty(vars.lambda)
	% interpolation at given dataCoords
	[controlPoints,~] = cp_interp(pb);
	% reparametrization
	t = reParametrize(t,vars.dataCoords);
else
	% fitting at given dataCoords
	[controlPoints,~] = cp_approx(pb,vars.lambda);
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

% evaluation at t of the curve and back-mapping to M
for i = 1:p
    tStart = tic;
    
    ti = t(i);
    
    sgmt = floor(ti)+1;
    if sgmt == n ; yTs = controlPoints(:,:,end); % last point
    elseif sgmt == ti+1 ; yTs = controlPoints(:,:,3*(sgmt-1)+1); % first or inner points
    else
		range = [3*(sgmt-1)+1:3*sgmt+1];
		
		% De casteljau on the two tangent spaces
		ti  = ti - sgmt + 1;
		yTs = euclideanBezier(controlPoints(:,:,range),ti);
    end
    
    y(i,:,:) = M.exp(rootPoint,yTs);
    tSpan(i) = toc(tStart);
end

if vars.display;
	fprintf('Blended spline evaluated in: '); 
	fprintf('%17.2f [s] \n\n',tSpan);
end

% Output tSpan if requested.
nout = max(nargout,1) - 1;
if nout == 2
	varargout{1} = curveVelocity(M,y,t);
  varargout{2} = tSpan;
elseif nout == 1
  varargout{1} = curveVelocity(M,y,t);
end

end


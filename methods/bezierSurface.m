function [y, varargout] = bezierSurface(varargin);
% Returns the value Z of the Bezier surface at coordinate X,Y.
%    The Bezier surface fits a set of data points DATAPOINTS
%    at coordinates DATACOORDS with a fitting focus LAMBDA.
%
%    Z = bezierSurface(M,DATAPOINTS,DATACOORDS,X,Y) returns Z, the 
%    value at coordinate [X,Y] of the Bezier surface interpolating the 
%    DATAPOINTS on a manifold M, associated with coordinates DATACOORDS
%    organized in a regular grid.
%    DATAPOINTS are stored in a matrix-cell of MxN entries and 
%    DATACOORDS are stored in a [MxNx2] matrix, where each line 
%    DATACOORDS(i,j,:) corresponds to the coordinate of the associated 
%    datapoint.
%       Another possibility is to provide data points in a vector-cell
%       of length m AND data coordinates in a [m x 2] matrix, 
%       representing a regular discretized grid.
% 
%    Z = bezierSurface(M,DATAPOINTS,[],X,Y) tries to interpolate the 
%    datapoints DATAPOINTS without the knowledge of the data coordinates.
%    In that case, the algorithm tries to associate data coordinates to 
%    the data points, on a regular grid. If not possible, an error is 
%    returned.
%
%    Z = bezierSurface(M,DATAPOINTS,DATACOORDS,X,Y,'options',OPTIONS)
%    permits to add some options to the system. The options are 
%    given hereunder.
%
%    [Z,TSPAN] = bezierSurface(M,DATAPOINTS,DATACOORDS,X,Y), returns 
%    also the time needed to compute each value of the fitting curve.
%
%    [Y,TSPAN,CONTROLPOINTS] = bezierSurface(M,DATAPOINTS,DATACOORDS,X,Y) 
%    returns also the control points of the surface.
%
%    Inputs: 
%      M is a manifold structure from Manopt - www.manopt.org.
%      DATAPOINTS must be a 2D-cell of dimension [m,n], where
%        and (m,n) is the number of data points in the x- and 
%        y-directions.
%      X, Y can be scalar values or a mesh grid of size p x q.
%      DATACOORDS is a matrix of size (m x n x 2), 
%        containing the coordinate at which each DATAPOINTS must be 
%        fitted. Default value : (i-1,j-1), i = 1:m, j = 1:n.
%      Optional arguments that can be:
%         'options' with a structure OPTIONS taking the following entries
%            options.curveGeneration to choose the algorithm to compute
%               the curve. Available possibilities are
%               * default | The default algorithm (recommended) is simple
%                         | and fast.
%               * tensor  | An algorithm that separates the generation
%                         | in x- and y-directions. Recommended when
%                         | m << n or n << m.
%               * optimal | A slow but more optimal way to compute the 
%                         | curve. Recommended when the points are far
%                         | from each other and when the exact optimal
%                         | transport is known on M.
%            options.curveReconstruction to choose the algorithm to
%            reconstruct the curve. Available possibilities are
%               * tensorx  | The default algorithm (recommended) is simple
%                          | and fast. It computes curves first in the
%                          | x-direction, then in the y-direction.
%               * tensory  | Same as tensorx, but inverts the order of
%                          | curve computation. Recommended as well.
%               * average  | Reconstruction algorithm based on the
%                          | karcher mean. Slower but sometimes more
%                          | accurate.
%              If you are extremely curious, here are deprecated versions 
%              that do not lead to c1 surfaces or even to non-continuous 
%              forms ;-).
%               * explog   | [Deprecated] Superfast method. Not continuous
%                          | on non-flat manifolds.
%               * karcher  | [Deprecated] Method based on the karcher
%                          | mean and the De Casteljau method.
%               * tensorx0 | Non-c1 version of tensorx (but c0).
%               * tensory0 | Non-c1 version of tensory (but c0).
%               * average0 | Non-c1 version of average (but c0).
%               
%    Outputs:
%      Y is a matrix of size [P,n_rows,n_cols] containing the values of the 
%        blended Bezier spline evaluated at the (P) times contained in T.
%      TSPAN (optional) is a value of the time needed to compute the surface.
%      CONTROLPOINTS (optional) returns a cell filled with the control points
%        of the Bezier surface.
%
% Original author: 
%   Pierre-Yves Gousenbourger, Apr. 12, 2019
% Contributors: 
%
% Change log:
% 	Apr. 12, 2019 (PYG):
%      First version
%   Apr. 16, 2019 (PYG):
%      Added the options for reconstruction and defense code
%      Added output of control points.

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
addRequired(ip,'dataCoords');
addRequired(ip,'t1')
addRequired(ip,'t2')
addOptional(ip,'display',false); % no display by default
  defOptions.curveReconstruction = 'tensorx';
  defOptions.curveGeneration = 'default';
addOptional(ip,'options',[]); % no default entry.

% parse the parser as a structure.
parse(ip, varargin{:});
vars = ip.Results;

% defense code
% dataPoints are given in the right display; 
if isvector(vars.dataPoints)
  assert(isRegularGrid(vars.dataCoords),'Data points and data coords are stored in an unstructured way. I can''t initialize the algorithm. Please provide at least a regular grid of coordinates.');
  assert(length(vars.dataPoints) == size(vars.dataCoords,1),'dataCoords and dataPoints must be of compatible size. Check the doc.');
  lengthX = length(unique(vars.dataCoords(:,1)));
  lengthY = length(unique(vars.dataCoords(:,2)));
  % sort the coordinates
  [~,indX] = sort(vars.dataCoords(:,1));
  vars.dataCoords = vars.dataCoords(indX,:);
  vars.dataPoints = vars.dataPoints(indX);
  [~,indY] = sort(vars.dataCoords(:,2));
  vars.dataCoords = vars.dataCoords(indY,:);
  vars.dataPoints = vars.dataPoints(indY);
  % create the matrices
  vars.dataPoints = reshape(vars.dataPoints,[lengthX,lengthY]);
  vars.dataCoords = reshape(vars.dataCoords,[lengthX,lengthY,2]);
end
assert(ismatrix(vars.dataPoints),'dataPoints must be stored in a 2d-cell.');
% t1 and t2 are the same size
assert(size(vars.t1,1) == size(vars.t2,1) && size(vars.t1,2) == size(vars.t2,2),'t1 and t2 must be the same size');
if ~isempty(vars.dataCoords)
  % dataCoords are size-compatible with dataPoints
  assert(size(vars.dataCoords,1) == size(vars.dataPoints,1) && size(vars.dataCoords,2) == size(vars.dataPoints,2), 'dataCoords and dataPoints must be of compatible size. Check the doc.');
  % t1 and t2 are not out of the dataCoords range
  assert(min(vars.t1(:)) >= min(min(vars.dataCoords(:,:,1))) && max(vars.t1(:)) <= max(max(vars.dataCoords(:,:,1))),'Please choose t1 inside the range of dataCoords(:,:,1).');
  assert(min(vars.t2(:)) >= min(min(vars.dataCoords(:,:,2))) && max(vars.t2(:)) <= max(max(vars.dataCoords(:,:,2))),'Please choose t2 inside the range of dataCoords(:,:,2).');
else
  % t1 and t2 are not out of default range
  maxT1 = size(vars.dataPoints,1)-1;
  maxT2 = size(vars.dataPoints,2)-1;
  assert(min(vars.t1(:)) >= 0 && max(vars.t1(:)) <= maxT1,['Please choose t1 between 0 and ',num2str(maxT1)]);
  assert(min(vars.t2(:)) >= 0 && max(vars.t2(:)) <= maxT2,['Please choose t2 between 0 and ',num2str(maxT2)]);
end

% fill in options
if ~isfield(vars.options,'curveReconstruction')
  vars.options.curveReconstruction = defOptions.curveReconstruction;
end
if ~isfield(vars.options,'curveGeneration')
  vars.options.curveGeneration = defOptions.curveGeneration;
end


% create the problem structure
pb.M = vars.M;
pb.space = pb.M.name();
pb.dataPoints = vars.dataPoints;
dimX       = size(vars.t1);
dimM       = size(vars.dataPoints{1,1});
t1 = vars.t1(:);
t2 = vars.t2(:);
p = size(t1);


% info display
if vars.display; 
	disp('=====================================================');
	fprintf('Bezier surface fitting tool: %26s\n',problem.M.name());
	disp('-----------------------------------------------------');
	tStart = tic;
	disp('Curve computation...');
end

% control points generation
switch vars.options.curveGeneration
case 'default'; generation = @control_points_simple_generation_2d;
case 'tensor';  generation = @control_points_double_tensorization;
case 'optimal'; generation = @control_points_generation;
otherwise; error('No such curveGeneration option value is known');
end
if ~isempty(vars.dataCoords)
    % interpolation at given dataCoords
    warning('The domain must be reparametrized. This is a transparent operation, no worries...');
    dcX = squeeze(vars.dataCoords(:,1,1));
    dcY = squeeze(vars.dataCoords(1,:,2));
    
    t1 = reParametrize(t1,dcX);
    t2 = reParametrize(t2,dcY);
end
pb = generation(pb);  


if vars.display;
	fprintf('Bezier surface computed in: ');  
	fprintf('%21.2f [s]\n\n',toc(tStart));
	disp('Curve evaluation...');
end

% Here, reconstruction method
switch vars.options.curveReconstruction
case 'tensorx';  curve = @(pb,X,Y) curve_reconstruction_double_bezier_c1(pb,X,Y,'direction','x-direction');
case 'tensory';  curve = @(pb,X,Y) curve_reconstruction_double_bezier_c1(pb,X,Y,'direction','y-direction');
case 'average';  curve = @curve_reconstruction_mean_all_c1;
case 'explog';   curve = @curve_reconstruction_de_casteljau_exp;
case 'karcher';  curve = @curve_reconstruction_de_casteljau_mean;
case 'tensorx0'; curve = @(pb,X,Y) curve_reconstruction_double_bezier(pb,X,Y,'direction','x-direction');
case 'tensory0'; curve = @(pb,X,Y) curve_reconstruction_double_bezier(pb,X,Y,'direction','y-direction');
case 'average0'; curve = @curve_reconstruction_mean_all;
otherwise; error('No such curveReconstruction option value is known');
end

tStart = tic;
pb = curve(pb,t1,t2);
y = squeeze(reshape(pb.curve,[dimX dimM]));

tSpan = toc(tStart);

if vars.display;
	fprintf('Bezier surface evaluated in: '); 
	fprintf('%17.2f [s] \n\n',tSpan);
end

% Output tSpan of cp if requested.
nout = max(nargout,1) - 1;
if nout >= 1
  varargout{1} = tSpan;
  if nout == 2
    varargout{2} = pb.controlPoints;
  else
    error('Too many output arguments');
  end
end

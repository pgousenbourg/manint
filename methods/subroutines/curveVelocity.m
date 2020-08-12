function [v, varargout] = curveVelocity(varargin);
% CURVEVELOCITY returns the value V of the velocity of the manifold-
%    valued curve GAMMA at a given time T. The velocity is given as a 
%    tangent vector at GAMMA(T). The velocity is computed as a forward
%    finite difference of the curve. The velocity at the last point of
%    GAMMA is the reverted finite difference between the last and 
%    before-last points.
%
%    V = CURVEVELOCITY(M,GAMMA,T) returns V, the value of the 
%    velocity of GAMMA at time T, on a manifold M. T are the values of
%    the parameter at which GAMMA is evaluated.
%
%    [V,TSPAN] = CURVEVELOCITY(M,DATAPOINTS,T) returns also the time 
%    needed to compute each value of the velocity.
%
%    Inputs: 
%      M is a manifold structure from Manopt - www.manopt.org.
%      GAMMA must be a matrix [n_rows,n_cols,n], where (n_rows,n_cols) 
%        is the dimension of the search space, and n is the number of 
%        points in GAMMA.
%      T is a vector of length n, containing the coordinate at which 
%        each value of GAMMA is evaluated.
%
%    Outputs:
%      V is a matrix of size [n_rows,n_cols,p] containing the values of 
%        the velocity of the curve GAMMA, at T. Note that V(:,:,i)
%        belongs to the Tangent Space of M at GAMMA(:,:,i).
%      TSPAN (optional) is a vector that contains the time spent to compute
%        the fitting curve at times contained in T.
%
% Original author: 
%   Pierre-Yves Gousenbourger, Nov. 26, 2018
% Contributors: 
%
% Change log:
% 	Nov. 26, 2018 (PYG):
%      First version

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
addRequired(ip,'M');
addRequired(ip,'gamma');
addRequired(ip,'t');

% parse the parser as a structure.
parse(ip, varargin{:});
vars = ip.Results;

% defense code
assert(length(size(vars.gamma)) == 3, 'GAMMA must be a 3-matrix. See help for more details.');
assert(size(vars.gamma,3) == length(vars.t),'Evaluation times and points on GAMMA must be the same.');

% variables
t     = vars.t;
gamma = vars.gamma;
[n_rows,n_cols,n] = size(gamma);

% output
v     = zeros(n_rows,n_cols,n);
tSpan = zeros(1,n);

% evaluation at t of the curve
for i = 1:n-1
    tStart   = tic;
    h        = abs(t(i+1) - t(i));
    v(:,:,i) = vars.M.log(gamma(:,:,i),gamma(:,:,i+1))./h;
    tSpan(i) = toc(tStart);
end

% evaluation of the end velocity
tStart   = tic;
h        = abs(t(end) - t(end-1));
v(:,:,end) = - vars.M.log(gamma(:,:,end),gamma(:,:,end-1))./h;
tSpan(end) = toc(tStart);

% Output tSpan if requested.
nout = max(nargout,1) - 1;
if nout == 1
   varargout{1} = tSpan;
end

end

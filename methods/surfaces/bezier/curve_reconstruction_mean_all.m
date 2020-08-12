function pb = curve_reconstruction_mean_all(varargin)
% Reconstructs the Bezier curve based on a karcher mean.
%
% function pb = curve_reconstruction_mean_all(pb,X,Y)
%    Reconstrucs the Bezier curve based on the control points stored in
%    the structure pb, at times (t1,t2). The piecewize Bezier surface is
%    piecewise cubic.
%
% Input: pb, the structure of the interpolating problem, containing
% 		   pb.dataPoints    (2-cell with data points)
% 		   pb.degree        (the degree of the curve: default 3)
%        pb.controlPoints (the control points).
%        X, the meshgrid of times in the X direction (pxq).
%        Y, the meshdrid of times in the Y direction (pxq).
%           
% Output: pb, the structure of the approximation problem, edited
%         with the field pb.curve in which the data are stored in a
%         (p,q,m,n) matrix.
%
% Original author: 
% 	Pierre-Yves Gousenbourger, Mar. 30, 2015.
% Contributors: 
%	  Paul Striewski, Jul. 14, 2015.
% Change log:
% 	Nov. 07, 2018 (PYG) - Integration of the Manopt framework.

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
  addRequired(ip,'pb');
  addRequired(ip,'X');
  addRequired(ip,'Y');
  addOptional(ip,'direction','x-direction'); % x-direction by default

  % parse the parser as a structure.
  parse(ip, varargin{:});
  vars = ip.Results;
  
  % variables
  pb = vars.pb;
  X  = vars.X;
  Y  = vars.Y;
  direction = vars.direction; 

  % parameters
  controlPoints = pb.controlPoints;
  [m,n]         = size(pb.dataPoints);
  [mCp,nCp]     = size(pb.controlPoints);
  if ~isfield(pb,'patchM'); pb.patchM = m-1; end
  if ~isfield(pb,'patchN'); pb.patchN = n-1; end
  if ~isfield(pb,'degree'); pb.degree = floor((mCp + 2)/m); end
  deg           = pb.degree;
  [r,c]         = size(pb.dataPoints{1,1});
  [mX,nX]       = size(X);
  M             = pb.M;
  
  
  % defense
  assert(mCp == 3*pb.patchM + 1,'Not enough control points');
  assert(nCp == 3*pb.patchN + 1,'Not enough control points');
  assert(sum(size(X) == size(Y))==2,'Matrices X and Y must be of the same dimension');
  assert((max(max(X)) <= pb.patchM) && (min(min(X) >= 0)),'times out of range');
  assert((max(max(Y)) <= pb.patchN) && (min(min(Y) >= 0)),'times out of range');
  
  % preallocation
  curve = zeros(mX,nX,r,c);
    
	
	% Waitbar just because it is fun ^^.
	h = waitbar(0,'Initialization of the reconstruction algoritm...');
	total_steps = mX*nX;
	perc = 0;
		
	% For each value of (X,Y)
	for k = 1:mX
	for l = 1:nX
		% find the patch
		patchX = max(1,ceil(X(k,l))) - 1;
		patchY = max(1,ceil(Y(k,l))) - 1;
		
		% reparametrize between 0 and 1
		t1 = X(k,l) - patchX;
		t2 = Y(k,l) - patchY;
		
		% control points
		rangeX = patchX*deg + 1 : (patchX+1)*deg + 1;
		rangeY = patchY*deg + 1 : (patchY+1)*deg + 1;
    
    % reduction
		cp = controlPoints(rangeX,rangeY);
		
		% Solution
    y = decasteljau2d_mean(M,cp,t1,t2);
		curve(k,l,:,:) = y;
		
		% Waitbar
		perc = perc + 1/total_steps;
		waitbar(perc,h,sprintf(...
			['Processing De Casteljau algorithm. Complete: %.f %%.'],...
			perc*100));
	end
	end
  
  close(h);

  %% Store
  pb.curve = curve;
end

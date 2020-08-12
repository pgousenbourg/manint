function b = bki(M,dataPoints, A, k, i, threshold)
% Computes the control point b_k^i.
%
% function b = bki(M, dataPoints, A, k, i, threshold)
%    returns b, a single control points based on A, in the segment i
%    and position k. The coefficients of A are cut at a given threshold 
%    around the diagonal.
%
% Input: M          (manifold structure of Manopt)
%        dataPoints (1-cell with data points)
%        A          (matrix defining the linear combinations of the
%                    dataPoints used for computing b_k^i)
%        i          (the segment)
%        k          (the position within the segment, between 0 and 3)
%        threshold  (the threshold applied in the lincomb).
%        
% Original author: 
% 	Pierre-Yves Gousenbourger, Nov. 04, 2015.
% Contributors: 
%	Paul Striewski.
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
	
	% Parameters
	m 	  = length(dataPoints);
	op 	  = [0,0,1,1];
	Tdata = cell(size(dataPoints));	% Projected dataPoints
	
	% Reference point.
	pref = dataPoints{i+op(k+1)};		% Closest point
	xbound = max(1,i-threshold-1):min(m,i+threshold+2);
	
	% Projection on tangent space
	for o = xbound
		Tdata{o} = M.log(pref,dataPoints{o}); 
	end
	
	alpha = cell(2,1);
	alpha{1} = compute_alpha(Tdata, A, i, threshold);
	alpha{2} = compute_alpha(Tdata, A, i+1, threshold);
	
	bk = 	((3-k)/3)	.*alpha{1} ...
		+ 	(k/3)		.*alpha{2};
		
	b 	= M.exp(pref, bk);
end

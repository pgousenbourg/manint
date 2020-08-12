function b = bklij(M, dataPoints, Am, An, k, l, i, j, threshold)
% Computes the control point b_kl^ij.
%
% function b = bklij(M, dataPoints, Am, An, k, l, i, j, threshold)
%    returns b, a single control points based on Am and An, in the patch
%    (i,j), and position (k,l) within the patch. The coefficients
%    of Am and An are cut at a given threshold around the diagonal.
%
% Input: M          (manifold structure of Manopt)
%        dataPoints (2-cell with data points)
%        Am, An     (matrices defining the linear combinations of the
%                    dataPoints used for computing b_kl^ij)
%        i,j        (the patch)
%        k,l        (the position within the patch, between 0 and 3)
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
	[m,n] = size(dataPoints);
	op 	 = [0,0,1,1];
	Tdata 	 = cell(size(dataPoints));	% Projected p
	
	% Reference point and projection to its tangent space.
	pref = dataPoints{i+op(k+1),j+op(l+1)};
	xbound = max(1,i-threshold-1):min(m,i+threshold+2);
	ybound = max(1,j-threshold-1):min(n,j+threshold+2);
	
	
	for o = xbound
	for p = ybound
		Tdata{o,p} = M.log(pref,dataPoints{o,p}); 
	end
	end
	
	% Step 1:
	% -------
	% compute [beta_(m,n-d:n+d)]
	% compute [beta_(m,n-d+1:n+d+1)]
	% compute [beta_(m,n-d+1:n+d+1)]
	% compute [beta_(m,n-d+1:n+d+1)]
	%
	% Step 2:
	% ------- 
	% use this to compute
	% a(m,n), a(m,n+1), a(m+1,n), a(m+1,n+1)
	beta = cell(n,1);
	alpha = cell(2,2);
	for s = 0:1
	for r = 0:1
		bound = max(1,j-threshold+r-1):min(n,j+threshold+r+1);
		% Step 1
		for o = bound
			beta{o} = compute_alpha(Tdata(:,o), Am, i+s, threshold);
		end
		% Step 2
		alpha{s+1,r+1} = compute_alpha(beta, An, j+r, threshold);
	end
	end
	
	bkl = 	(((3-k)/3)*((3-l)/3))	.*alpha{1,1} ...
		+ 	(((3-k)/3)*(l/3))		.*alpha{1,2} ...
		+ 	((k/3)*((3-l)/3))		.*alpha{2,1} ...
		+ 	((k/3)*(l/3))			.*alpha{2,2} ;
		
	b 	= M.exp(pref, bkl);
end

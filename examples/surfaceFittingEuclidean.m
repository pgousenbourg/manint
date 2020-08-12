% surfaceFittingEuclidean.m
% 
% This file is an example of surface fitting on the Euclidean space.
% It shows the basic use of the functions provided in the toolbox of
% my thesis, ie.,
%    * blendedSurface
%    * bezierSurface
%    * localSurface
%    * geoSurface
% 
% Author: Pierre-Yves Gousenbourger.
% Version: May. 13, 2020
% log: May. 13, 2020 - PYG
%       First version

close all;
clear all;
clc;


% First of all, you must add the methods to the path.
% The toolbox is based on MANOPT. The simpliest way to install manopt is
% to add the complete folder.
% The methods provided in my thesis are in the folder /methods.
addpath(genpath([pwd,'/../manopt']));
addpath(genpath([pwd,'/../methods']));


% Second, one needs to specify several data points.
% For instance, one might want to approximate the following function
f = @(x,y) cos(2.*x) + sin(3.*y.^2) + cos(2.*x + 5.*y);
% based on a small number of data points.

% We thus generate a set of coordinates ordered on a regular grid (this
% is not necessary for blended surfaces, but well for bezier surface and
% piecewise geodesic surfaces).
[x,y] = meshgrid(linspace(0,1,4),linspace(0,1,6));
x = x(:); y = y(:);
% and store the data points in a cell of size sD
dataCoords = [x,y];
for i = 1:length(x)
  dataPoints{i,1} = [x(i),y(i),f(x(i),y(i))];
end
clear x y;

% Third, we must prepare the algorithms by specifying the manifold and
% eventually rootPoints (for the blending method, for instance).
dim = [1,3];
M   = euclideanfactory(dim);
% We decide to have three patches in the x-direction, and two in y,
% situated on f.
[x,y] = meshgrid(linspace(0,1,4),linspace(0,1,3));
x = x(:); y = y(:);
rootCoords = [x,y];
for i = 1:length(x)
  rootPoints{i,1} = [x(i),y(i),f(x(i),y(i))];
end
clear x y;


% Finally, we can set the coordinates at which we want to reconstruct f:
[X,Y] = meshgrid(linspace(0,1,41),linspace(0,1,31));
[m,n] = size(X);


% The algorithms are performed here.
% There are four different algorithms for surface fitting and interpolation:
% The blendedSurface (fitting), the bezierSurface (interpolation), the 
% localSurface (fitting on only one tangent space) and geoSurface (geodesic
% interpolation).
% 
% The algorithms my have several options. We refer to the documentation
% included in the help of each method for more details.
options.rootPoints = rootPoints;
options.rootCoords = rootCoords;
lambda = NaN; % interpolation

% Here, we compute the different surface fittings.
% Interpolation with blended surfaces...
blendedSurf = blendedSurface(M,dataPoints,dataCoords,X,Y,'lambda',lambda,'options',options);
% Interpolation with bezier surfaces...
bezierSurf  = bezierSurface(M,dataPoints,dataCoords,X,Y,'options',options);
% Interpolation with local surfaces in the tangent space of the first root point...
localSurf   = localSurface(M,dataPoints,dataCoords,X,Y,'rootPoint',rootPoints{1});
% Interpolation with geodesic surfaces
geoSurf     = geoSurface(M,dataPoints,dataCoords,X,Y);



% After some post-processing, we can plot the results.
blendedSurf = reshape(blendedSurf,m,n,3);
bezierSurf  = reshape(bezierSurf,m,n,3);
localSurf   = reshape(localSurf,m,n,3);
geoSurf     = reshape(geoSurf,m,n,3);

dP = cell2mat(dataPoints);

% Plot of the blended surface
figure;
surf(blendedSurf(:,:,1),blendedSurf(:,:,2),blendedSurf(:,:,3)); hold on;
plot3(dP(:,1),dP(:,2),dP(:,3),'.r','MarkerSize',30);

% Plot of the Bezier surface
figure;
surf(bezierSurf(:,:,1),bezierSurf(:,:,2),bezierSurf(:,:,3)); hold on;
plot3(dP(:,1),dP(:,2),dP(:,3),'.r','MarkerSize',30);

% Plot of the local surface
figure;
surf(localSurf(:,:,1),localSurf(:,:,2),localSurf(:,:,3)); hold on;
plot3(dP(:,1),dP(:,2),dP(:,3),'.r','MarkerSize',30);

% Plot of the geo surface
figure;
surf(geoSurf(:,:,1),geoSurf(:,:,2),geoSurf(:,:,3)); hold on;
plot3(dP(:,1),dP(:,2),dP(:,3),'.r','MarkerSize',30);

% surfaceFittingSphere.m
% 
% This file is an example of surface fitting on the sphere S(2).
% It shows the basic use of the functions provided in the toolbox of
% my thesis. Here we show only an example of use of blendedSurface.
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
% Let select a range of values of phi and theta:
t1 = linspace(0,1,41)';
t2 = linspace(0,1,31)';
delta = linspace(pi/3,-pi/3,41);
theta = linspace(0,2*pi/3,31) + pi;
% and a selection of indexes.
indX   = [17 15 31 32 7 20 18 26 29 30 11 27 26 6 4 20 39 13 23 9];
indY   = [10 25 18 17 28 8 23 23 11 17 2 1 16 24 28 4 17 14 31 10];
% The mapping function to the sphere is.
f = @(xx,yy) [cos(xx).*cos(yy), cos(xx).*sin(yy), sin(xx)];


% We generate a set of coordinates as well as the data points:
dataCoords = [t1(indX),t2(indY)];
for i = 1:length(indX)
  dataPoints{i,1} = f(delta(indX(i))',theta(indY(i))');
end


% Third, we must prepare the algorithms by specifying the manifold and
% eventually rootPoints (for the blending method, for instance).
M   = spherefactory(3);
% We decide to have three patches in the x-direction, and two in y,
% situated on f.
[x,y] = meshgrid(linspace(0,1,4),linspace(0,1,3));
[delta,theta] = meshgrid(linspace(pi/3,-pi/3,4),linspace(pi,5*pi/3,3));
x = x(:); y = y(:);
delta = delta(:); theta = theta(:);
rootCoords = [x,y];
for i = 1:length(x)
  rootPoints{i,1} = f(delta(i),theta(i));
end



% Finally, we can set the coordinates at which we want to reconstruct f:
X = linspace(0,1,100);
Y = 0.5*(sin(4*pi*X)+1);
m = length(X);

% The algorithm is performed here.
options.rootPoints = rootPoints;
options.rootCoords = rootCoords;
lambda = 100; % fitting

blendedSurf = blendedSurface(M,dataPoints,dataCoords,X,Y,'lambda',lambda,'options',options);


% After some post-processing, one can plot the result on a sphere
bS = blendedSurf;
dP = cell2mat(dataPoints);
rP = cell2mat(rootPoints);

figure;
% sphere
[XS,YS,ZS] = sphere(20);
surf(XS*0.99,YS*0.99,ZS*0.99,'FaceColor',[255 215 0]/255,'FaceAlpha',0.5,'EdgeAlpha',0.5); 
hold on; axis equal off; view(-40,5);
% blended surface
plot3(bS(:,1),bS(:,2),bS(:,3),'lineWidth',2); hold on;
% points
plot3(dP(:,1),dP(:,2),dP(:,3),'.r','MarkerSize',30);
plot3(rP(:,1),rP(:,2),rP(:,3),'.g','MarkerSize',30);

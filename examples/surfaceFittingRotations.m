% surfaceFittingRotations.m
% 
% This file is an example of surface fitting on the group SO(3).
% It shows the basic use of the functions provided in the toolbox of
% my thesis. Here we show only an example of use of blendedSurface,
% bezierSurface, geoSurface and localSurface, for points on a regular
% grid.
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
addpath(genpath([pwd,'/helpers'])); % some tools for data generation and plots


% Second, one needs to specify several data points and data coordinates.
% The data points will here also be the root points. We choose to have
% 3 patches in x and in y
patches = 3;
sampling = 2;
[x,y] = meshgrid(0:patches,0:patches);
x = x(:); y = y(:);
dataCoords = [x,y];
for i = 1:length(x)
  dataPoints{i} = SO3matrix(); % random data points
end
clear x y;

% Third, we must prepare the algorithms by specifying the manifold and
% the coordinates at which the surface must be evaluated.
M   = rotationsfactory(3);
[X,Y] = meshgrid(linspace(0,patches,patches*sampling+1),linspace(0,patches,patches*sampling+1));

% The algorithm is performed here.
options.rootPoints = dataPoints;
options.rootCoords = dataCoords;

% No lambda is provided : interpolation
S{1,1} = blendedSurface(M,dataPoints,dataCoords,X,Y,'options',options);
S{2,1} = bezierSurface(M,dataPoints,dataCoords,X,Y);
S{3,1} = localSurface(M,dataPoints,dataCoords,X,Y);
S{4,1} = geoSurface(M,dataPoints,dataCoords,X,Y);


% The plot procedure is a bit long, but the result is nice ;-)
dimM = [3,3];
dimX = size(X);
toy = model3d('helpers/3Dobjects/Toy car N130810.3DS');
D = 450;
fplot = @(rot, type, place) draw_so3_object(rot,toy,type,place);
titles = {'blended','bezier','local','geodesics'};
for k = 1:size(S,1)

  Z = S{k};
  figure;
  % Surface
  for i = 1:dimX(1)
  for j = 1:dimX(2)
      fplot(reshape(Z(i,j,:,:),dimM),'normal',[D*(i-1) , D*(j-1), 0]);
  end
  end

  % Data points
  data = reshape(dataPoints,patches+1,patches+1);
  for i = 1:4;
  for j = 1:4;
    fplot(data{i,j},'interp',[D*sampling*(i-1), D*sampling*(j-1), 0]);
  end
  end
  light; lighting phong;
  title(titles{k});
end

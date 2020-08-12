% curveFittingRotations.m
% 
% This file is an example of curve fitting on the group SO(3).
% It shows the basic use of the functions provided in the toolbox of
% my thesis. Here we show only a basic example of use of blendedCurve,
% bezierCurve, TScubic and interpGeo. More options are possible,
% be we refer to the help of each function for that.
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
% It is possible to not specify any coordinates: just leave the data
% coordinates empty and the algorithm will consider by default that 
% data point i is at coordinate t=i.
sampling = 3;
nData = 6;
dataCoords = linspace(0,1,nData);
for i = 1:nData
  dataPoints{i,1} = SO3matrix(); % random data points
end

% Third, we must prepare the algorithms by specifying the manifold and
% the coordinates at which the surface must be evaluated.
M = rotationsfactory(3);
t = linspace(0,1,(nData-1)*sampling+1);
lambda = 10; % fitting

% The fitting curve is produced here (four different methods)
S{1,1} = blendedCurve(M,dataPoints,dataCoords,t,'lambda',lambda);
S{2,1} = bezierCurve(M,dataPoints,dataCoords,t,'lambda',lambda);
S{3,1} = TScubic(M,dataPoints,dataCoords,t,'lambda',lambda); % by default the rootPoint is in the middle; you can optionally choose it.
S{4,1} = interpGeo(M,dataPoints,dataCoords,t); % interpolation only


% The plot procedure is a bit long, but the result is worth it ;-)
dimM = [3,3];
dimX = length(t);
toy = model3d('helpers/3Dobjects/Toy car N130810.3DS');
D = 450;
fplot = @(rot, type, place) draw_so3_object(rot,toy,type,place);
titles = {'blended','bezier','local','geodesics'};

for k = 1:size(S,1)
  Z = S{k};
  figure;
  % Surface
  for i = 1:dimX
      fplot(reshape(Z(i,:,:),dimM),'normal',[D*(i-1),0,0]);
  end

  % Data points
  for i = 1:length(dataPoints);
    fplot(dataPoints{i},'interp',[D*sampling*(i-1),0,0]);
  end
  light; lighting phong;
  title(titles{k});
end

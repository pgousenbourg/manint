% curveFittingSphere.m
% 
% This file is an example of curve fitting on the Sphere S(2).
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


% Second, one needs to specify several data points and data coordinates.
% It is possible to not specify any coordinates: just leave the data
% coordinates empty and the algorithm will consider by default that 
% data point i is at coordinate t=i.
sampling = 3;
f = @(theta,phi) [cos(theta).*sin(phi) , sin(theta).*sin(phi), cos(phi)];
phi     = [pi/3 ; 4*pi/7 ; 4*pi/7 ; pi/3] - 2*pi/12; 
theta   = [pi/2 ; pi/4 ; 3*pi/4 ; pi/2] + 10*pi/12;
dataPoints = mat2cell(f(theta,phi),ones(1,length(phi)),3);
nData   = length(dataPoints);
dataCoords = linspace(0,1,nData);

% Third, we must prepare the algorithms by specifying the manifold and
% the coordinates at which the surface must be evaluated.
M = spherefactory(3);
t = linspace(0,1,(nData-1)*30);
lambda = 10; % fitting

% The fitting curve is produced here (four different methods)
S{1,1} = blendedCurve(M,dataPoints,dataCoords,t,'lambda',lambda);
S{2,1} = bezierCurve(M,dataPoints,dataCoords,t,'lambda',lambda);
S{3,1} = TScubic(M,dataPoints,dataCoords,t,'lambda',lambda); % by default the rootPoint is in the middle; you can optionally choose it.
S{4,1} = interpGeo(M,dataPoints,dataCoords,t); % interpolation only


% Now we just plot the results on the same sphere
titles = {'blended','bezier','local','geodesic'};
dP = cell2mat(dataPoints);


figure;
% sphere
  [XS,YS,ZS] = sphere(20);
  surf(XS*0.99,YS*0.99,ZS*0.99,'FaceColor',[255 215 0]/255,'FaceAlpha',0.5,'EdgeAlpha',0.5,'HandleVisibility','off'); 
  hold on; axis equal off; view(-40,5);
% curves
  for i = 1:4
    Z = reshape(S{i},[length(t),3]);
    plot3(Z(:,1),Z(:,2),Z(:,3),'lineWidth',2); hold on;
  end
% data points
plot3(dP(:,1),dP(:,2),dP(:,3),'.r','MarkerSize',30);
legend(titles);

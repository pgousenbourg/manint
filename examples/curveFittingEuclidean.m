% curveFittingEuclidean.m
% 
% This file is an example of curve fitting on the Euclidean space R2.
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
nData = 6;
dataCoords = linspace(0,1,nData);
dataPoints = [0 1 1 2 2 3 ; 0 0 1 1 2 2]';
dataPoints = mat2cell(dataPoints,ones(1,nData),2);


% Third, we must prepare the algorithms by specifying the manifold and
% the coordinates at which the surface must be evaluated.
M = euclideanfactory([1,2]);
t = linspace(0,1,100);
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
% curves
for i = 1:4
  subplot(2,2,i);
  Z = reshape(S{i},[length(t),2]);
  plot(Z(:,1),Z(:,2),'lineWidth',2); hold on;
  plot(dP(:,1),dP(:,2),'.r','MarkerSize',30);
  title(titles{i});
end

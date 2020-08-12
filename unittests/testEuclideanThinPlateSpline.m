% Unit tests: Euclidean Thin Plate Spline.
% 
% Files tested:
%   * thinPlateSplineGeneration.m
%   * thinPlateSplineReconstruction.m
% 
% This file is part of the project "bezierfitting" with B. Wirth from
% uni-muenster and PY. Gousenbourger from UCLouvain.
% 
% Author: Pierre-Yves Gousenbourger.
% Version: Dec. 17, 2019
% log: Dec. 17, 2019 - PYG
%       First version

close all;

addpath(genpath([pwd,'/../methods']));


disp('Unit test on Euclidean Thin Plate Spline');
success  = 0;
failures = 0;


% ======================================================================
% Tests 1 to 8: a flat solution. 
% Tests of the structure of the inputs and outputs.
% ======================================================================
mdp = 5;
ndp = 4;

% data
const = 1;
for i = 1:mdp*ndp
  dataPoints{i} = const;
end
for j = 1:ndp
for i = 1:mdp
  dataCoords(i+(j-1)*mdp,:) = [i-1,j-1];
end
end
[X,Y] = meshgrid(linspace(0,mdp-1,(mdp-1)*10+1),linspace(0,ndp-1,(ndp-1)*10+1));
  
% --- test: constant solution, interpolation.
fprintf('  test constant solution...');

  repr = thinPlateSplineGeneration(dataPoints,dataCoords);
  tps   = thinPlateSplineReconstruction(X,Y,repr,dataCoords);
  
  % Conditions
  passed = (norm(tps - const) == 0);
  if passed
    fprintf(' passed!\n');
    success = success + 1;
  else
    fprintf(2,' error!\n');
    failures = failures + 1;
  end
  
  clear passed tps;
  

% --- test: thinPlateSplineGeneration with wrong dataPoints formatting
fprintf('  detection of wrong dataPoints formatting...');

try
  repr = thinPlateSplineGeneration(cell2mat(dataPoints),dataCoords);
  % This part should never be reached. If so, then it means that the
  % error wasn't thrown...
  fprintf(2,' error!\n');
  failures = failures + 1;
catch ME
  if strcmp(ME.identifier,'surfaceFitting:dataNotInCell')
    fprintf(' passed!\n');
    success = success + 1;
  else
    fprintf(2,' error!\n');
    failures = failures + 1;
  end
end

clear tps;

% --- test: thinPlateSplineGeneration with wrong dataCoords dimension
fprintf('  detection of wrong dataCoords dimension...');

passed = 1;

try
  repr = thinPlateSplineGeneration(dataPoints,dataCoords(:,1));
  % This part should never be reached. If so, then it means that the
  % error wasn't thrown...
  passed = 0;
catch ME
  if strcmp(ME.identifier,'surfaceFitting:coordsNotIn2D')
    passed = passed*1;
  else
    passed = 0;
  end
end
try
  repr = thinPlateSplineGeneration(dataPoints,dataCoords');
  % This part should never be reached. If so, then it means that the
  % error wasn't thrown...
  fprintf(2,' error!\n');
  failures = failures + 1;
catch ME
  if strcmp(ME.identifier,'surfaceFitting:coordsNotIn2D')
    passed = passed*1;
  else
    passed = 0;
  end
end

if passed
  fprintf(' passed!\n');
  success = success + 1;
else
  fprintf(2,' error!\n');
  failures = failures + 1;
end

clear passed;

% --- test: thinPlateSplineGeneration with incompatible dataPoints and
% dataCoords
fprintf('  detection of incompatible dataPoints/dataCoords...');

try
  repr = thinPlateSplineGeneration(dataPoints,dataCoords(1:end-1,:));
  % This part should never be reached. If so, then it means that the
  % error wasn't thrown...
  fprintf(2,' error!\n');
  failures = failures + 1;
catch ME
  if strcmp(ME.identifier,'surfaceFitting:dimCheck')
    fprintf(' passed!\n');
    success = success + 1;
  else
    fprintf(2,' error!\n');
    failures = failures + 1;
  end
end

clear tps

% --- test: thinPlateSplineReconstruction with wrong X-Y inputs
fprintf('  detection of X and Y of different sizes...');

try
  tps   = thinPlateSplineReconstruction(X,1,repr,dataCoords);
  % This part should never be reached. If so, then it means that the
  % error wasn't thrown...
  fprintf(2,' error!\n');
  failures = failures + 1;
catch ME
  if strcmp(ME.identifier,'surfaceFitting:dimCheck')
    fprintf(' passed!\n');
    success = success + 1;
  else
    fprintf(2,' error!\n');
    failures = failures + 1;
  end
end

clear tps;



% ======================================================================
% Tests: affine solution (interpolation)
% ======================================================================

mdp = 5;
ndp = 4;

% data
for j = 1:ndp
for i = 1:mdp
  dataCoords(i+(j-1)*mdp,:) = [i-1,j-1];
  dataPoints{i+(j-1)*mdp} = i+j;
end
end
[X,Y] = meshgrid(linspace(0,mdp-1,(mdp-1)*10+1),linspace(0,ndp-1,(ndp-1)*10+1));

% --- test: affine solution, interpolation.
fprintf('  test affine symmetric solution...');

  repr = thinPlateSplineGeneration(dataPoints,dataCoords);
  tps   = thinPlateSplineReconstruction(X,Y,repr,dataCoords);
  
  passed = 1;
  for i = 1:mdp
  for j = 1:ndp
    passed = passed*(tps((j-1)*10+1,(i-1)*10+1) == i+j);
  end
  end
  % Conditions
  %passed = (norm(tps - const) == 0);
  if passed
    fprintf(' passed!\n');
    success = success + 1;
  else
    fprintf(2,' error!\n');
    failures = failures + 1;
  end
  
  clear passed;
  clear mdp ndp const dataPoints dataCoords X Y a d tps;

% ======================================================================
% Tests: nonlinear solutions (interpolation)
% ======================================================================

mdp = 5;
ndp = 4;

% data
for j = 1:ndp
for i = 1:mdp
  dataCoords(i+(j-1)*mdp,:) = [i-1,j-1];
  val = cos(2*i) + sin(3*j^2) + cos(2*i+5*j);
  dataPoints{i+(j-1)*mdp} = val;
end
clear val;
end
[X,Y] = meshgrid(linspace(0,mdp-1,(mdp-1)*10+1),linspace(0,ndp-1,(ndp-1)*10+1));


% --- test: quadratic solution, interpolation.
fprintf('  test nonlinear solution...');

  repr = thinPlateSplineGeneration(dataPoints,dataCoords);
  tps   = thinPlateSplineReconstruction(X,Y,repr,dataCoords);
  
  % conditions
  passed = 1;
  tol = 1e-12;
  for i = 1:mdp
  for j = 1:ndp
    val = cos(2*i) + sin(3*j^2) + cos(2*i+5*j);
    passed = passed*(tps((j-1)*10+1,(i-1)*10+1) - val < tol);
  end
  end
  if passed
    fprintf(' passed!\n');
    success = success + 1;
  else
    fprintf(2,' error!\n');
    failures = failures + 1;
  end

  clear val tol
  clear mdp ndp const dataPoints dataCoords X Y a d tps;


% ======================================================================
% Tests: nonregular grid nonlinear solutions
% ======================================================================

m = 5;
n = 4;

% data
[X,Y] = meshgrid(linspace(0,m-1,(m-1)*10+1),linspace(0,n-1,(n-1)*10+1));

xSelec = [17 15 31 32 7 20 18 26 29 30 11 27 26 6 4 20 39 13 23 9];
ySelec = [10 25 18 17 28 8 23 23 11 17 2 1 16 24 28 4 17 14 31 10];
ndp = length(ySelec);
dataCoords = [X(1,xSelec)',Y(ySelec,1)];
for i = 1:ndp
  x = dataCoords(i,1);
  y = dataCoords(i,2);
  val = cos(2*x) + sin(3*y^2) + cos(2*x+5*y);
  dataPoints{i} = val;
end
clear val;

% --- test: non uniform grid
fprintf('  test nonuniform grid...');

  repr = thinPlateSplineGeneration(dataPoints,dataCoords);
  tps   = thinPlateSplineReconstruction(X,Y,repr,dataCoords);
  
  %figure;
  %dp = cell2mat(dataPoints);
  %plot3(dataCoords(:,1),dataCoords(:,2),dp,'.r','MarkerSize',30);
  %hold on;
  %surf(X,Y,tps);
  
  % conditions
  passed = 1;
  tol = 1e-12;
  for i = 1:ndp
    x = dataCoords(i,1);
    y = dataCoords(i,2);
    val = cos(2*x) + sin(3*y^2) + cos(2*x+5*y);
    
    passed = passed*(tps(ySelec(i),xSelec(i)) - val < tol);
  end
  if passed
    fprintf(' passed!\n');
    success = success + 1;
  else
    fprintf(2,' error!\n');
    failures = failures + 1;
  end

  clear val tol
  clear m n ndp xSelect ySelec dataPoints dataCoords X Y a d tps;


% ======================================================================
% Tests: fitting tests (regular grid, constant solution)
% ======================================================================

mdp = 5;
ndp = 4;
lambda = 0;

% data
const = 1;
for i = 1:mdp*ndp
  dataPoints{i} = const;
end
for j = 1:ndp
for i = 1:mdp
  dataCoords(i+(j-1)*mdp,:) = [i-1,j-1];
end
end

% --- test on constant solution with lambda to 0
fprintf('  test constant solution with lambda=0...');

  repr = thinPlateSplineGeneration(dataPoints,dataCoords,lambda);
  a = repr.a;
  
  % Conditions
  passed = (a(1,:,:) == 1 && a(2,:,:) == 0 && a(3,:,:) == 0);
  if passed
    fprintf(' passed!\n');
    success = success + 1;
  else
    fprintf(2,' error!\n');
    failures = failures + 1;
  end
  
  clear passed;
  clear mdp ndp const lambda dataPoints dataCoords a d;

% ======================================================================
% Tests: fitting tests (regular grid, nonlinear solution, lambda = 0)
% ======================================================================

mdp = 5;
ndp = 4;
lambda = 0;

% data
for j = 1:ndp
for i = 1:mdp
  dataCoords(i+(j-1)*mdp,:) = [i-1,j-1];
  val = cos(2*i) + sin(3*j^2) + cos(2*i+5*j);
  dataPoints{i+(j-1)*mdp} = val;
end
clear val;
end
[X,Y] = meshgrid(linspace(0,mdp-1,(mdp-1)*10+1),linspace(0,ndp-1,(ndp-1)*10+1));


% --- test: quadratic solution, interpolation.
fprintf('  test nonlinear solution with lambda=0...');

  repr = thinPlateSplineGeneration(dataPoints,dataCoords,lambda);
  tps   = thinPlateSplineReconstruction(X,Y,repr,dataCoords);
  
  % conditions
  passed = 1;
  tol = 1e-12;
  for i = 1:mdp
  for j = 1:ndp
    val = cos(2*i) + sin(3*j^2) + cos(2*i+5*j);
    passed = passed*(tps((j-1)*10+1,(i-1)*10+1) - val < tol);
  end
  end
  if passed
    fprintf(' passed!\n');
    success = success + 1;
  else
    fprintf(2,' error!\n');
    failures = failures + 1;
  end

  clear passed;
  clear lambda val a d tps tol;

% --- test: quadratic solution, interpolation.
fprintf('  test nonlinear solution with lambda=1e-2...');

lambda = 1e-2;

try
  repr = thinPlateSplineGeneration(dataPoints,dataCoords,lambda);
  tps   = thinPlateSplineReconstruction(X,Y,repr,dataCoords);
  
  %figure;
  %plot3(dataCoords(:,1),dataCoords(:,2),cell2mat(dataPoints),'.r','MarkerSize',30);
  %hold on;
  %surf(X,Y,tps);
  
  fprintf(' passed!\n');
  success = success + 1;
catch  
  fprintf(2,' error!\n');
  failures = failures + 1;
end

  clear val tol
  clear mdp ndp const dataPoints dataCoords X Y a d tps lambda;


% ======================================================================
% Tests: test of matrix dataPoint (nonlinear, interpolation)
% ======================================================================

mdp = 5;
ndp = 4;

% data
for j = 1:ndp
for i = 1:mdp
  dataCoords(i+(j-1)*mdp,:) = [i-1,j-1];
  val = cos(2*i) + sin(3*j^2) + cos(2*i+5*j);
  dataPoints{i+(j-1)*mdp} = [i-1,j-1,val];
end
dataPoints = dataPoints(:);
clear val;
end
[X,Y] = meshgrid(linspace(0,mdp-1,(mdp-1)*10+1),linspace(0,ndp-1,(ndp-1)*10+1));


% --- test: interpolation on vector spaces.
fprintf('  test nonlinear solution on vector space...');

  repr = thinPlateSplineGeneration(dataPoints,dataCoords);
  tps   = thinPlateSplineReconstruction(X,Y,repr,dataCoords);
  
  dp = cell2mat(dataPoints);
  Z  = reshape(tps,[size(X),3]);
  
  %figure
  %plot3(dp(:,1),dp(:,2),dp(:,3),'.r','MarkerSize',30);
  %hold on;
  %surf(Z(:,:,1),Z(:,:,2),Z(:,:,3));
  
  % conditions
  passed = 1;
  tol = 1e-12;
  for i = 1:mdp
  for j = 1:ndp
    val = cos(2*i) + sin(3*j^2) + cos(2*i+5*j);
    passed = passed*(tps((j-1)*10+1,(i-1)*10+1,1) - i < tol);
    passed = passed*(tps((j-1)*10+1,(i-1)*10+1,2) - j < tol);
    passed = passed*(tps((j-1)*10+1,(i-1)*10+1,3) - val < tol);
  end
  end
  if passed
    fprintf(' passed!\n');
    success = success + 1;
  else
    fprintf(2,' error!\n');
    failures = failures + 1;
  end
  clear val tol passed

% --- test: interpolation on matrix spaces.
fprintf('  test nonlinear solution on matrix space...');

% data
for j = 1:ndp
for i = 1:mdp
  dataCoords(i+(j-1)*mdp,:) = [i-1,j-1];
  val = cos(2*i) + sin(3*j^2) + cos(2*i+5*j);
  dataPoints{i+(j-1)*mdp} = repmat([i-1,j-1,val],[3,1]);
end
clear val;
end
dataPoints = dataPoints(:);
  
  repr = thinPlateSplineGeneration(dataPoints,dataCoords);
  tps   = thinPlateSplineReconstruction(X,Y,repr,dataCoords);
  
  dp = cell2mat(dataPoints);
  Z  = reshape(tps(:,:,1,:),[size(X),3]);
  
  %figure;
  %plot3(dp(:,1),dp(:,2),dp(:,3),'.r','MarkerSize',30);
  %hold on;
  %surf(Z(:,:,1),Z(:,:,2),Z(:,:,3));
  
  % conditions
  passed = 1;
  tol = 1e-12;
  for i = 1:mdp
  for j = 1:ndp
    val = cos(2*i) + sin(3*j^2) + cos(2*i+5*j);
    passed = passed*(tps((j-1)*10+1,(i-1)*10+1,1,1) - i < tol);
    passed = passed*(tps((j-1)*10+1,(i-1)*10+1,2,1) - i < tol);
    passed = passed*(tps((j-1)*10+1,(i-1)*10+1,3,1) - i < tol);
    passed = passed*(tps((j-1)*10+1,(i-1)*10+1,1,2) - j < tol);
    passed = passed*(tps((j-1)*10+1,(i-1)*10+1,2,2) - j < tol);
    passed = passed*(tps((j-1)*10+1,(i-1)*10+1,3,2) - j < tol);
    passed = passed*(tps((j-1)*10+1,(i-1)*10+1,1,3) - val < tol);
    passed = passed*(tps((j-1)*10+1,(i-1)*10+1,2,3) - val < tol);
    passed = passed*(tps((j-1)*10+1,(i-1)*10+1,3,3) - val < tol);
  end
  end
  if passed
    fprintf(' passed!\n');
    success = success + 1;
  else
    fprintf(2,' error!\n');
    failures = failures + 1;
  end
  
  clear val tol passed
  clear mdp ndp const dataPoints dataCoords a d dp Z repr;

% ======================================================================

nbTests = success + failures;
fprintf('\n  Number of tests:%i\n',nbTests);
fprintf('   -- %i success(es) (%1.0f %%)\n',success,100*success/nbTests);
fprintf('   -- %i failure(s) (%1.0f %%)\n\n',failures,100*failures/nbTests);

% Unit tests: liftings
% 
% Files tested:
%   * pointsLifting.m
%   * basicLifting.m (to change)
%   * coarseLifting.m (todo)
%   * transportedLifting.m (todo)
% 
% This file is part of the project "bezierfitting" with B. Wirth from
% uni-muenster and PY. Gousenbourger from UCLouvain.
% 
% Author: Pierre-Yves Gousenbourger.
% Version: Dec. 19, 2019
% log: Dec. 19, 2019 - PYG
%       First version

close all;

addpath(genpath([pwd,'/../methods']));
addpath(genpath([pwd,'/../manopt']));

disp('Unit tests on the lifting methods');
success  = 0;
failures = 0;


% ======================================================================
% Points Lifting
% ======================================================================

% --- basics test: one point to lift
fprintf('  test POINTSLIFTING: basicTests 1...');
  
  M    = euclideanfactory(1); %
  x    = 1;
  y{1} = 1;
  
  passed = (cell2mat(pointsLifting(M,x,y)) == 0);
  if passed; fprintf(' passed!\n'); success = success + 1; 
  else; fprintf(2,' error!\n'); failures = failures + 1;
  end

  clear passed x y;

% --- basics test: multiple points to lift (vector)
fprintf('  test POINTSLIFTING: basicTests 2 (vector)...');
  
  M    = euclideanfactory(1); %
  x    = 1;
  y{1} = 1; y{2} = 2;
  result = cell2mat(pointsLifting(M,x,y));
  
  passed = (sum(result - [0 1])==0);
  if passed; fprintf(' passed!\n'); success = success + 1; 
  else; fprintf(2,' error!\n'); failures = failures + 1;
  end

  clear passed x y result;
  
% --- basics test: multiple points to lift (matrix)
fprintf('  test POINTSLIFTING: basicTests 3 (matrix)...');
  
  M    = euclideanfactory(1); %
  x    = 1;
  y{1,1} = 11; y{1,2} = 12; y{2,1} = 21; y{2,2} = 22;
  result = cell2mat(pointsLifting(M,x,y));
  
  passed = (sum(result(:) - [10;20;11;21]) == 0);
  if passed; fprintf(' passed!\n'); success = success + 1; 
  else; fprintf(2,' error!\n'); failures = failures + 1;
  end

  clear passed x y result;
  
% --- basics test: multiple points to lift (multidimensional matrix)
fprintf('  test POINTSLIFTING: basicTests 4 (tensors)...');
  
  M    = euclideanfactory(1); %
  x    = 1;
  y{1,1,1} = 111;
  y{2,1,1} = 211;
  y{1,2,1} = 121;
  y{2,2,1} = 221;
  y{1,1,2} = 112;
  y{2,1,2} = 212;
  y{1,2,2} = 122;
  y{2,2,2} = 222;
  expected = cell2mat(y(:))-1;
  result = cell2mat(pointsLifting(M,x,y));
  
  passed = (sum(result(:) - expected) == 0);
  if passed; fprintf(' passed!\n'); success = success + 1; 
  else; fprintf(2,' error!\n'); failures = failures + 1;
  end

  clear passed x y result expected;
  
% --- test assertions
fprintf('  test POINTSLIFTING: rootsNotInVector...');
M    = euclideanfactory(1); %
passed = 1;
try % too many entries in weights
  pointsLifting(M,{1,1;2,2},1);
  passed = 0;
catch ME
  if ~strcmp(ME.identifier,'surfaceFitting:vectorCell'); passed = 0; end
end

  if passed; fprintf(' passed!\n'); success = success + 1; 
  else; fprintf(2,' error!\n'); failures = failures + 1;
  end

clear passed M;

% ======================================================================

nbTests = success + failures;
fprintf('\n  Number of tests:%i\n',nbTests);
fprintf('   -- %i success(es) (%1.0f %%)\n',success,100*success/nbTests);
fprintf('   -- %i failure(s) (%1.0f %%)\n\n',failures,100*failures/nbTests);

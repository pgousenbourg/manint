% Unit tests: blendSurfaces_old
% 
% Files tested:
%   * blendSurfWeights.m
%   * blendSurfaces.m
% 
% This file is part of the project "bezierfitting" with B. Wirth from
% uni-muenster and PY. Gousenbourger from UCLouvain.
% 
% Author: Pierre-Yves Gousenbourger.
% Version: Dec. 18, 2019
% log: Dec. 18, 2019 - PYG
%       First version

close all;

addpath(genpath([pwd,'/../methods']));
addpath(genpath([pwd,'/../manopt']));

disp('Unit tests on the averaging methods');
success  = 0;
failures = 0;



% ======================================================================
% blendSurfWeights
% ======================================================================

% --- test basics
fprintf('  test SURFACEWEIGHTS helper: basicTest 1...');
  
  X = 0; Y = 0;
  rootCoords = [0 0; 1 0; 0 1; 1 1];
  expectedWeights = [1 0 0 0];
  
  if (sum(blendSurfWeights(X,Y,rootCoords) - expectedWeights) == 0); 
    fprintf(' passed!\n'); 
    success = success + 1; 
  else 
    fprintf(2,' error!\n'); 
    failures = failures + 1;
  end

  clear X Y rootCoords expectedWeights;


% --- test basics
fprintf('  test SURFACEWEIGHTS helper: basicTest 2...');
  
  t = [0 0; 0 1; 1 0; 1 1];
  X = t(:,1); Y = t(:,2);
  rootCoords = [0 0; 1 0; 0 1; 1 1];
  expectedWeights = [1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1];
  
  if (sum(sum(cell2mat(blendSurfWeights(X,Y,rootCoords)) - expectedWeights)) == 0); 
    fprintf(' passed!\n'); 
    success = success + 1; 
  else 
    fprintf(2,' error!\n'); 
    failures = failures + 1;
  end

  clear t X Y rootCoords expectedWeights;


% --- test: no rootCoords given
fprintf('  test SURFACEWEIGHTS helper: no rootCoords...');

  t = [0 0; 0 1; 1 0; 1 1];
  X = t(:,1); Y = t(:,2);
  expectedWeights = [1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1];
  
  if (sum(sum(cell2mat(blendSurfWeights(X,Y)) - expectedWeights)) == 0); 
    fprintf(' passed!\n'); 
    success = success + 1; 
  else 
    fprintf(2,' error!\n'); 
    failures = failures + 1;
  end

  clear t X Y expectedWeights;


% --- test: non-unitary rootCoords
fprintf('  test SURFACEWEIGHTS helper: non-unitary rootCoords...');
  
  t = [0 0; 0 1; 1 0; 1 1]./2;
  X = t(:,1); Y = t(:,2);
  rootCoords = [0 0; 1 0; 0 1; 1 1]./2;
  expectedWeights = [1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1];
  
  if (sum(sum(cell2mat(blendSurfWeights(X,Y,rootCoords)) - expectedWeights)) == 0); 
    fprintf(' passed!\n'); 
    success = success + 1; 
  else 
    fprintf(2,' error!\n'); 
    failures = failures + 1;
  end

  clear t X Y rootCoords expectedWeights;


% --- test: non-unitary rootCoords (2)
fprintf('  test SURFACEWEIGHTS helper: scaling...');

  t = [1 1]./4;
  X = t(:,1); Y = t(:,2);
  rootCoords = [0 0; 1 0; 0 1; 1 1]./2;
  expectedWeights = 0.25*ones(1,4);
  
  if (sum(sum(blendSurfWeights(X,Y,rootCoords) - expectedWeights)) == 0); 
    fprintf(' passed!\n'); 
    success = success + 1; 
  else 
    fprintf(2,' error!\n'); 
    failures = failures + 1;
  end

  clear t X Y rootCoords expectedWeights;


% --- input test: rootCoords flipped by 90°
fprintf('  test SURFACEWEIGHTS helper: flippedBasis...');

t = [0 0; 0 1; 1 0; 1 1];
X = t(:,1); Y = t(:,2);
rootCoords = [0 0; 0 1; 1 0; 1 1]; % here the grid is flipped
try
  blendSurfWeights(X,Y,rootCoords);
  % This part should never be reached.
  passed = 0;
catch ME
  if strcmp(ME.identifier,'surfaceWeights:gridMisAlignment'); passed = 1;
  else passed = 0;
  end
end

if passed; fprintf(' passed!\n'); success = success + 1; 
else fprintf(2,' error!\n'); failures = failures + 1;
end

clear passed t X Y rootCoords

% --- input test: X and Y of different sizes
fprintf('  test SURFACEWEIGHTS helper: dimCheck(1)...');

X = 0;
Y = [0 0 1];
try
  blendSurfWeights(X,Y);
  % This part should never be reached.
  passed = 0;
catch ME
  if strcmp(ME.identifier,'surfaceFitting:dimCheck'); passed = 1;
  else passed = 0;
  end
end

if passed; fprintf(' passed!\n'); success = success + 1; 
else fprintf(2,' error!\n'); failures = failures + 1;
end

clear passed X Y

% --- input test: rootPoints wronly sized
fprintf('  test SURFACEWEIGHTS helper: dimCheck(2)...');
X = 0;
Y = 0; 
rootPoints = [1 1];
try
  blendSurfWeights(X,Y,rootPoints);
  % This part should never be reached.
  passed = 0;
catch ME
  if strcmp(ME.identifier,'surfaceFitting:dimCheck'); passed = 1;
  else passed = 0;
  end
end

if passed; fprintf(' passed!\n'); success = success + 1; 
else fprintf(2,' error!\n'); failures = failures + 1;
end

clear rootPoints X Y passed;

% --- advanced test: weights sum to one
fprintf('  test SURFACEWEIGHTS helper: w sum to 1...');

t = [linspace(0,1,5);linspace(0,1,5)]';
X = t(:,1); Y = t(:,2);
rootCoords = [0 0; 1 0; 0 1; 1 1];
w = cell2mat(blendSurfWeights(X,Y,rootCoords));
passed = 1;
passed = passed*(sum(sum(w,2) - 1) < 1e-12);

t = [linspace(0,1,4);linspace(0,1,4)]';
X = t(:,1); Y = t(:,2);
rootCoords = [0 0; 1 0; 0 1; 1 1];
w = cell2mat(blendSurfWeights(X,Y,rootCoords));
%passed = 1;
passed = passed*(sum(sum(w,2) - 1) < 1e-12);

if passed; fprintf(' passed!\n'); success = success + 1; 
else fprintf(2,' error!\n'); failures = failures + 1;
end

clear passed t X Y rootCoords ;


% ==========================================
% BlendSurfaces_test: On the Euclidean space
% ==========================================

% --- basics test: just one point (R^1)
fprintf('  test blendSurface: basicTests 1...');
  
  %
  M = euclideanfactory(1);
  rootPoints{1}    = M.rand(); 
  rootPoints = repmat(rootPoints,1,4);
  
  tangentValues{1} = 2; 
  tangentValues = repmat(tangentValues,1,4);
  
  weights  = [1 0 0 0];
  
  expected = rootPoints{1}+2;
  result   = blendSurfaces(M,rootPoints,tangentValues,weights);
  
  passed = (result == expected);
  if passed; fprintf(' passed!\n'); success = success + 1; 
  else; fprintf(2,' error!\n'); failures = failures + 1;
  end

  clear passed rootPoints tangentValues weights expected result;
  
% --- basics test: just one point (R^2)
fprintf('  test blendSurface: basicTests 2...');
  
  %
  M = euclideanfactory(1,2);
  rootPoints{1} = M.rand(); 
  rootPoints    = repmat(rootPoints,1,4);
  
  tangentValues{1} = [2,2]; 
  tangentValues    = repmat(tangentValues,1,4);
  
  weights  = [1 0 0 0];
  
  expected = rootPoints{1}+2;
  result   = blendSurfaces(M,rootPoints,tangentValues,weights);
  
  passed = (sum(result - expected) == 0);
  
  if passed; fprintf(' passed!\n'); success = success + 1; 
  else; fprintf(2,' error!\n'); failures = failures + 1;
  end

  clear passed rootPoints tangentValues weights expected result;

% --- basics test: just one point (R^(3 x 2))
fprintf('  test blendSurface: basicTests 3...');
  
  M = euclideanfactory(3,2);
  rootPoints{1} = M.rand(); 
  rootPoints    = repmat(rootPoints,1,4);
  
  tangentValues{1} = repmat(2,3,2); 
  tangentValues    = repmat(tangentValues,1,4);
  
  weights  = [1 0 0 0];
  
  expected = rootPoints{1}+2;
  result   = blendSurfaces(M,rootPoints,tangentValues,weights);
  
  passed = (sum(result - expected) == 0);
  
  if passed; fprintf(' passed!\n'); success = success + 1; 
  else; fprintf(2,' error!\n'); failures = failures + 1;
  end

  clear passed rootPoints tangentValues weights expected result;
  
  
% --- basics test: multiple points (same dimension)
fprintf('  test blendSurface: multipoints...');
  
  M = euclideanfactory(1);
  rootPoints = cell(1,4);
  for i = 1:4
    rootPoints{i} = i;
  end
  
  tangentValues{1} = [1 1 1 1]';
  tangentValues = repmat(tangentValues,1,4);
  
  weights{1} = [1 1 1 1]';
  weights{2} = [0 0 0 0]';
  weights{3} = weights{2};
  weights{4} = weights{2};
  
  expected = 2*ones(4,1);
  result   = blendSurfaces(M,rootPoints,tangentValues,weights);
  
  passed = (sum(expected - result) == 0);
  
  if passed; fprintf(' passed!\n'); success = success + 1; 
  else; fprintf(2,' error!\n'); failures = failures + 1;
  end

  clear M rootPoints tangentValues weights i passed expected result;


% --- basics test: multiple points (only one line of weights)
fprintf('  test blendSurface: multipoints (one weight)...');

  M = euclideanfactory(1);
  rootPoints = cell(1,4);
  for i = 1:4
    rootPoints{i} = i;
  end
  
  tangentValues{1} = [1 1 1 1]';
  tangentValues = repmat(tangentValues,1,4);
  
  weights  = [1 0 0 0];
  
  expected = 2*ones(4,1);
  result   = blendSurfaces(M,rootPoints,tangentValues,weights);
  
  passed = (sum(expected - result) == 0);

  if passed; fprintf(' passed!\n'); success = success + 1; 
  else; fprintf(2,' error!\n'); failures = failures + 1;
  end

  clear M rootPoints tangentValues weights i passed expected result;
  
% --- basics test: optional handle "average"
fprintf('  test blendSurface: handle average...');
  
  average = @(N,r,w) tensorMean(N,r,w,'orientation',2);
  
  M = euclideanfactory(1);
  rootPoints = cell(1,4);
  for i = 1:4
    rootPoints{i} = i;
  end
  
  tangentValues{1} = [1 1 1 1]';
  tangentValues = repmat(tangentValues,1,4);
  
  weights  = [1 0 0 0];
  
  expected = 2*ones(4,1);
  result   = blendSurfaces(M,rootPoints,tangentValues,weights,average);
  
  passed = (sum(expected - result) == 0);
  
  if passed; fprintf(' passed!\n'); success = success + 1; 
  else; fprintf(2,' error!\n'); failures = failures + 1;
  end

  clear M rootPoints tangentValues weights i passed expected result;



% --- test assertions
M = euclideanfactory(1);

fprintf('  test blendSurface: Roots Not In Cell...');

passed = 1;
try
  blendSurfaces(M,1,cell(1),1);
  passed = 0;
catch ME
  if ~strcmp(ME.identifier,'surfaceFitting:dataNotInCell'); passed = 0; end
end
if passed; fprintf(' passed!\n'); success = success + 1; 
else; fprintf(2,' error!\n'); failures = failures + 1;
end

clear passed;


% --- test assertions
fprintf('  test blendSurface: TangentValues Not In Cell...');

passed = 1;
try
  blendSurfaces(M,cell(1),1,1);
  passed = 0;
catch ME
  if ~strcmp(ME.identifier,'surfaceFitting:dataNotInCell'); passed = 0; end
end
if passed; fprintf(' passed!\n'); success = success + 1; 
else; fprintf(2,' error!\n'); failures = failures + 1;
end

clear passed;


% --- test assertions
fprintf('  test blendSurface: length of entries...');

passed = 1;
try % no 4 entries in tangentValues
  blendSurfaces(M,cell(1),cell(1),1);
  passed = 0;
catch ME
  if ~strcmp(ME.identifier,'blendSurfaces:fourPoints'); passed = 0; end
end
try % no 4 entries in rootPoints
  blendSurfaces(M,cell(1,4),cell(1),1);
  passed = 0;
catch ME
  if ~strcmp(ME.identifier,'blendSurfaces:fourPoints'); passed = 0; end
end
try % no 4 entries in weights
  blendSurfaces(M,cell(1,4),cell(1,4),cell(1));
  passed = 0;
catch ME
  if ~strcmp(ME.identifier,'blendSurfaces:fourPoints'); passed = 0; end
end
if passed; fprintf(' passed!\n'); success = success + 1; 
else; fprintf(2,' error!\n'); failures = failures + 1;
end

clear passed;


% --- test assertions
fprintf('  test blendSurface: tests on weights (cell of vectors)...');

passed = 1;
rootPoints = cell(1,4);

tangentValues{1} = zeros(5,1);
tangentValues = repmat(tangentValues,1,4);

try % weights incompatible with TV
  weights{1} = zeros(2,1);
  weights = repmat(weights,1,4);

  blendSurfaces(M,rootPoints,tangentValues,weights);
  passed = 0;
catch ME
  if ~strcmp(ME.identifier,'surfaceFitting:dimCheck'); passed = 0; end
end
if passed; fprintf(' passed!\n'); success = success + 1; 
else; fprintf(2,' error!\n'); failures = failures + 1;
end

clear passed rootPoints tangentValues weights;

% --- test assertions
fprintf('  test blendSurface: tests on weights (cell of matrices)...');

passed = 1;
rootPoints = cell(1,4);

tangentValues{1} = zeros(5,2,4);
tangentValues = repmat(tangentValues,1,4);

try % weights incompatible with TV
  weights{1} = zeros(2,5);
  weights = repmat(weights,1,4);

  blendSurfaces(M,rootPoints,tangentValues,weights);
  passed = 0;
catch ME
  if ~strcmp(ME.identifier,'surfaceFitting:dimCheck'); passed = 0; end
  clear weights;
end
try % weights incompatible with TV
  weights{1} = zeros(5,3);
  weights = repmat(weights,1,4);

  blendSurfaces(M,rootPoints,tangentValues,weights);
  passed = 0;
catch ME
  if ~strcmp(ME.identifier,'surfaceFitting:dimCheck'); passed = 0; end
  clear weights;
end
if passed; fprintf(' passed!\n'); success = success + 1; 
else; fprintf(2,' error!\n'); failures = failures + 1;
end

clear passed rootPoints tangentValues weights;

% --- test assertions
fprintf('  test blendSurface: tests on weights (vector)...');

passed = 1;
try % not enough weights
  blendSurfaces(M,cell(1,4),cell(1,4),[1,2]);
  passed = 0;
catch ME
  if ~strcmp(ME.identifier,'blendSurfaces:fourPoints'); passed = 0; end
end
if passed; fprintf(' passed!\n'); success = success + 1; 
else; fprintf(2,' error!\n'); failures = failures + 1;
end

clear passed;


% --- test assertions
fprintf('  test blendSurface: tests on weights (?!$&)...');

passed = 1;
try % not enough weights
  blendSurfaces(M,cell(1,4),cell(1,4),[1,2;3,4]);
  passed = 0;
catch
end
if passed; fprintf(' passed!\n'); success = success + 1; 
else; fprintf(2,' error!\n'); failures = failures + 1;
end

clear passed;

% ======================================================================

nbTests = success + failures;
fprintf('\n  Number of tests:%i\n',nbTests);
fprintf('   -- %i success(es) (%1.0f %%)\n',success,100*success/nbTests);
fprintf('   -- %i failure(s) (%1.0f %%)\n\n',failures,100*failures/nbTests);

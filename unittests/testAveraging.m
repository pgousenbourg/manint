% Unit tests: Methods for Averaging
% 
% Files tested:
%   * tensorMean.m
%   * karcherMean.m
% 
% This file is part of the project "bezierfitting" with B. Wirth from
% uni-muenster and PY. Gousenbourger from UCLouvain.
% 
% Author: Pierre-Yves Gousenbourger.
% Version: Dec. 18, 2019
% log: Dec. 18, 2019 - PYG
%        First version
%      Jan. 06, 2020 - PYG
%        Integration of the new framework.
%      Jan. 07, 2020 - PYG
%        Tests for karcherMean.m

close all;

addpath(genpath([pwd,'/../methods']));
addpath(genpath([pwd,'/../manopt']));

disp('Unit tests on the averaging methods');
success  = 0;
failures = 0;


% ======================================================================
% tensorMean
% ======================================================================
% --- basics test: just one patch, one point to evaluate.
fprintf('  test tensorMean: basicTests 1...');
  
  M = euclideanfactory(3,2);
  
  points = cell(1,4);
  for i = 1:4; points{i} = i*ones(3,2); end
  weights = blendSurfWeights(0,0);
  expectedPoint = points{1};
  average = tensorMean(M,points,weights);
  if sum(sum(expectedPoint - average)) == 0; fprintf(' passed!\n'); success = success + 1; 
  else; fprintf(2,' error!\n'); failures = failures + 1;
  end

  clear points weights expectedPoint average;

% --- basics test: just one patch, one point, no weights
fprintf('  test tensorMean: basicTests 2...');
  
  points = cell(1,4);
  for i = 1:4; points{i} = i*ones(3,2); end
  expectedPoint = 2.5*ones(3,2);
  average = tensorMean(M,points);
  if sum(sum(expectedPoint - average)) == 0; fprintf(' passed!\n'); success = success + 1; 
  else; fprintf(2,' error!\n'); failures = failures + 1;
  end

  clear points expectedPoint average;

% --- basics test: just one patch, multiple points to evaluate.
fprintf('  test tensorMean: basicTests 3...');
  
  points = cell(1,4);
  for i = 1:4;
    for j = 1:10
      temp(j,:,:) = 10*(i-1)*ones(3,2) + j;
    end
    points{i} = temp;
  end
  clear temp;
  
  weights = cell2mat(blendSurfWeights(zeros(10,1),zeros(10,1)));
  expectedPoints = points{1};
  average = tensorMean(M,points,weights);
  if sum(sum(expectedPoints - average)) == 0; fprintf(' passed!\n'); success = success + 1; 
  else; fprintf(2,' error!\n'); failures = failures + 1;
  end

  clear points weights expectedPoints average;

% --- basics test: just one patch, multiple points, no weight
fprintf('  test tensorMean: basicTests 4...');
  
  points = cell(1,4);
  for i = 1:4;
    for j = 1:10
      temp(j,:,:) = 10*(i-1)*ones(3,2) + j;
    end
    points{i} = temp;
  end
  clear temp;
  
  expectedPoint = .25*(points{1} + points{2} + points{3} + points{4});
  average = tensorMean(M,points);
  if sum(sum(expectedPoint - average)) == 0; fprintf(' passed!\n'); success = success + 1; 
  else; fprintf(2,' error!\n'); failures = failures + 1;
  end

  clear points expectedPoint average;

% --- basic test: orientation changed, with weights
fprintf('  test tensorMean: basicTests 5...');
  
  points = cell(1,4);
  for i = 1:4;
    for j = 1:10
      temp(j,:,:) = 10*(i-1)*ones(3,2) + j;
    end
    points{i} = temp;
  end
  clear temp;
  weights = cell2mat(blendSurfWeights(zeros(10,1),zeros(10,1)));
  
  expectedPoint = points{1};
  average = tensorMean(M,points,weights,'orientation',2);
  if sum(sum(expectedPoint - average)) == 0; fprintf(' passed!\n'); success = success + 1; 
  else; fprintf(2,' error!\n'); failures = failures + 1;
  end

  clear points weights expectedPoint average;


% --- test assertions
fprintf('  test tensorMean: dataNotInCell...');

passed = 0;
try 
  tensorMean(M,1);
catch ME
  if strcmp(ME.identifier,'surfaceFitting:dataNotInCell'); passed = 1; end
end
if passed; fprintf(' passed!\n'); success = success + 1; 
else; fprintf(2,' error!\n'); failures = failures + 1;
end

clear passed;

% --- test assertions
fprintf('  test tensorMean: length of points...');

passed = 0;
try 
  tensorMean(M,cell(1,1));
catch ME
  if strcmp(ME.identifier,'blendSurf:fourEntries'); passed = 1; end
end
if passed; fprintf(' passed!\n'); success = success + 1; 
else; fprintf(2,' error!\n'); failures = failures + 1;
end

clear passed;


% --- test assertions
fprintf('  test tensorMean: size of points...');

passed = 0;
points{1} = ones(3,2);
points{2} = points{1};
points{3} = zeros(size(points{1}));
points{4} = 1;

try 
  tensorMean(M,points);
catch ME
  if strcmp(ME.identifier,'surfaceFitting:dimCheck'); passed = 1; end
end
if passed; fprintf(' passed!\n'); success = success + 1; 
else; fprintf(2,' error!\n'); failures = failures + 1;
end

clear passed points;



% --- test assertions
fprintf('  test tensorMean: dataNotInMatrix...');

passed = 0;
points{1} = ones(3,2);
points = repmat(points,1,4);
try 
  tensorMean(M,points,cell(1,1));
catch ME
  if strcmp(ME.identifier,'blendSurf:fourEntries'); passed = 1; end
end
if passed; fprintf(' passed!\n'); success = success + 1; 
else; fprintf(2,' error!\n'); failures = failures + 1;
end

clear passed points;

% --- test assertions
fprintf('  test tensorMean: length of weights...');

passed = 0;
points{1} = ones(3,2);
points = repmat(points,1,4);
try 
  tensorMean(M,points,2);
catch ME
  if strcmp(ME.identifier,'blendSurf:fourEntries'); passed = 1; end
end
if passed; fprintf(' passed!\n'); success = success + 1; 
else; fprintf(2,' error!\n'); failures = failures + 1;
end

clear passed points;


% --- test assertions

fprintf('  test tensorMean: dimCheck...');

passed = 0;
points{1} = ones(3,2);
points{2} = zeros(3,2);
points{3} = ones(3,2);
points{4} = zeros(3,2);

weights = repmat([1 1 1 1],2,1);

try 
  tensorMean(M,points,weights);;
catch ME
  if strcmp(ME.identifier,'surfaceFitting:dimCheck'); passed = 1; end
end
if passed; fprintf(' passed!\n'); success = success + 1; 
else; fprintf(2,' error!\n'); failures = failures + 1;
end

clear passed points weights


% --- test assertions

fprintf('  test tensorMean: orientationMismatch...');

points{1} = ones(3,2);
points{2} = zeros(3,2);
points{3} = ones(3,2);
points{4} = zeros(3,2);

passed = 1;
try 
  tensorMean(M,points,'orientation',3);
  passed = 0;
catch ME
  if ~strcmp(ME.identifier,'tensorMean:orientation'); passed = 0; end
end
try
  tensorMean(M,points,'orientation',0);
  passed = 0;
catch ME
  if ~strcmp(ME.identifier,'tensorMean:orientation'); passed = 0; end
end
if passed; fprintf(' passed!\n'); success = success + 1; 
else; fprintf(2,' error!\n'); failures = failures + 1;
end

clear passed points M;


% ======================================================================
% karcherMean
% ======================================================================
% --- basics test: just one patch, one point to evaluate.
fprintf('  test karcherMean: basicTests 1...');
  
  M = euclideanfactory(3,2);
  tol = 1e-10;
  
  points = cell(1,4);
  for i = 1:4; points{i} = i*ones(3,2); end
  weights = blendSurfWeights(0,0);
  expectedPoint = reshape(points{1},[3,2]);
  average = karcherMean(M,points,weights);
  if sum(sum(expectedPoint - average)) < tol; fprintf(' passed!\n'); success = success + 1; 
  else; fprintf(2,' error!\n'); failures = failures + 1;
  end

  clear points weights expectedPoint average M tol;

% --- basics test: just one patch, one point, no weights
fprintf('  test karcherMean: basicTests 2...');
  
  M = euclideanfactory(3,2);
  tol = 1e-10;
  
  points = cell(1,4);
  for i = 1:4; temp(1,:,:) = i*ones(3,2); points{i} = temp; end
  clear temp;
  expectedPoint = 2.5*ones(3,2);
  average = karcherMean(M,points);
  if sum(sum(expectedPoint - average)) < tol; fprintf(' passed!\n'); success = success + 1; 
  else; fprintf(2,' error!\n'); failures = failures + 1;
  end

  clear points expectedPoint average tol;

% --- basics test: just one patch, multiple points to evaluate.
fprintf('  test karcherMean: basicTests 3...');
  
  tol = 1e-10;
  points = cell(1,4);
  for i = 1:4;
    for j = 1:10
      temp(j,:,:) = 10*(i-1)*ones(3,2) + j;
    end
    points{i} = temp;
  end
  clear temp;
  
  weights = cell2mat(blendSurfWeights(zeros(10,1),zeros(10,1)));
  expectedPoints = points{1};
  average = karcherMean(M,points,weights);
  if sum(sum(sum(expectedPoints - average))) < tol; fprintf(' passed!\n'); success = success + 1; 
  else; fprintf(2,' error!\n'); failures = failures + 1;
  end

  clear points weights expectedPoints average tol;

% --- basics test: just one patch, multiple points, no weight
fprintf('  test karcherMean: basicTests 4...');
  
  tol = 1e-8;
  points = cell(1,4);
  for i = 1:4;
    for j = 1:10
      temp(j,:,:) = 10*(i-1)*ones(3,2) + j;
    end
    points{i} = temp;
  end
  clear temp;
  
  expectedPoint = .25*(points{1} + points{2} + points{3} + points{4});
  average = karcherMean(M,points);
  if sum(sum(sum(expectedPoint - average))) < tol; fprintf(' passed!\n'); success = success + 1; 
  else; fprintf(2,' error!\n'); failures = failures + 1;
  end

  clear points expectedPoint average tol;


% --- test assertions
fprintf('  test karcherMean: dataNotInCell...');

passed = 0;
try 
  karcherMean(M,1);
catch ME
  if strcmp(ME.identifier,'surfaceFitting:dataNotInCell'); passed = 1; end
end
if passed; fprintf(' passed!\n'); success = success + 1; 
else; fprintf(2,' error!\n'); failures = failures + 1;
end

clear passed;

% --- test assertions
fprintf('  test karcherMean: length of points...');

passed = 0;
try 
  karcherMean(M,cell(1,1));
catch ME
  if strcmp(ME.identifier,'blendSurf:fourEntries'); passed = 1; end
end
if passed; fprintf(' passed!\n'); success = success + 1; 
else; fprintf(2,' error!\n'); failures = failures + 1;
end

clear passed;


% --- test assertions
fprintf('  test karcherMean: size of points...');

passed = 0;
points{1} = ones(3,2);
points{2} = points{1};
points{3} = zeros(size(points{1}));
points{4} = 1;

try 
  karcherMean(M,points);
catch ME
  if strcmp(ME.identifier,'surfaceFitting:dimCheck'); passed = 1; end
end
if passed; fprintf(' passed!\n'); success = success + 1; 
else; fprintf(2,' error!\n'); failures = failures + 1;
end

clear passed points;



% --- test assertions
fprintf('  test karcherMean: dataNotInMatrix...');

passed = 0;
points{1} = ones(3,2);
points = repmat(points,1,4);
try 
  karcherMean(M,points,cell(1,1));
catch ME
  if strcmp(ME.identifier,'blendSurf:fourEntries'); passed = 1; end
end
if passed; fprintf(' passed!\n'); success = success + 1; 
else; fprintf(2,' error!\n'); failures = failures + 1;
end

clear passed points;

% --- test assertions
fprintf('  test karcherMean: length of weights...');

passed = 0;
points{1} = ones(3,2);
points = repmat(points,1,4);
try 
  karcherMean(M,points,2);
catch ME
  if strcmp(ME.identifier,'blendSurf:fourEntries'); passed = 1; end
end
if passed; fprintf(' passed!\n'); success = success + 1; 
else; fprintf(2,' error!\n'); failures = failures + 1;
end

clear passed points;


% --- test assertions

fprintf('  test karcherMean: dimCheck...');

passed = 0;
points{1} = ones(3,2);
points{2} = zeros(3,2);
points{3} = ones(3,2);
points{4} = zeros(3,2);

weights = repmat([1 1 1 1],2,1);

try 
  karcherMean(M,points,weights);;
catch ME
  if strcmp(ME.identifier,'surfaceFitting:dimCheck'); passed = 1; end
end
if passed; fprintf(' passed!\n'); success = success + 1; 
else; fprintf(2,' error!\n'); failures = failures + 1;
end

clear passed points weights


% ======================================================================

nbTests = success + failures;
fprintf('\n  Number of tests:%i\n',nbTests);
fprintf('   -- %i success(es) (%1.0f %%)\n',success,100*success/nbTests);
fprintf('   -- %i failure(s) (%1.0f %%)\n\n',failures,100*failures/nbTests);

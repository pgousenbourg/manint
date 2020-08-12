% Unit tests: Helpers
% 
% Files tested:
%   * rbf.m
%   * rootPos.m
%   * isRegularGrid.m
%   * sortCoords.m
%   * isParent.m
%   * isSemi.m
%   * findRoot.m
%   * computeLimits.m
%   * findNeighbors.m
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


disp('Unit tests on the helpers');
success  = 0;
failures = 0;


% ======================================================================
% RBF
% ======================================================================
fprintf('  test RBF helper: 0 is returned...');  
  if (rbf(0) == 0); fprintf(' passed!\n'); success = success + 1;
  else fprintf(2,' error!\n'); failures = failures + 1; end

fprintf('  test RBF helper: rbf(1)==0...');
  if (rbf(1) == 0); fprintf(' passed!\n'); success = success + 1;
  else fprintf(2,' error!\n'); failures = failures + 1; end

fprintf('  test RBF helper: rbf(1)==0...');
  if (rbf(2) == (1./(16*pi))*4*log(4)); 
    fprintf(' passed!\n'); success = success + 1;
  else fprintf(2,' error!\n'); failures = failures + 1; end


% ======================================================================
% rootPos
% ======================================================================
% --- test basics
fprintf('  test ROOTPOS helper: basicTests 1...');
  pos = [1 2 3 4];
  if (sum(rootPos(pos) - [1 1; 2 1; 1 2; 2 2]) == 0); 
    fprintf(' passed!\n'); 
    success = success + 1; 
  else 
    fprintf(2,' error!\n'); 
    failures = failures + 1;
  end

  clear pos;

% --- test basics
fprintf('  test ROOTPOS helper: basicTests 2...');
  i = [1 1 2 2];
  j = [1 2 1 2];
  if (sum(rootPos(i,j) - [1 3 2 4]') == 0); 
    fprintf(' passed!\n'); 
    success = success + 1; 
  else 
    fprintf(2,' error!\n'); 
    failures = failures + 1;
  end

  clear i j;

% --- test assertions

fprintf('  test ROOTPOS helper: input tests...');
passed = 1;

% i > 4
try
  rootPos(5);
  % This part should never be reached.
  passed = 0;
catch ME
  if ~strcmp(ME.identifier,'rootPos:inputType'); passed = 0; end
end

% i < 1
try
  rootPos(0);
  % This part should never be reached.
  passed = 0;
catch ME
  if ~strcmp(ME.identifier,'rootPos:inputType'); passed = 0; end
end

% double input, i < 1
try
  rootPos(0,1);
  % This part should never be reached.
  passed = 0;
catch ME
  if ~strcmp(ME.identifier,'rootPos:inputType'); passed = 0; end
end

% double input i > 2
try
  rootPos(3,1);
  % This part should never be reached.
  passed = 0;
catch ME
  if ~strcmp(ME.identifier,'rootPos:inputType'); passed = 0; end
end

% double input j < 1
try
  rootPos(1,0);
  % This part should never be reached.
  passed = 0;
catch ME
  if ~strcmp(ME.identifier,'rootPos:inputType'); passed = 0; end
end

% double input j > 2
try
  rootPos(1,3);
  % This part should never be reached.
  passed = 0;
catch ME
  if ~strcmp(ME.identifier,'rootPos:inputType'); passed = 0; end
end

if passed
  fprintf(' passed!\n'); 
  success = success + 1; 
else 
  fprintf(2,' error!\n'); 
  failures = failures + 1;
end

clear passed;


% ======================================================================
% isRegularGrid
% ======================================================================
% --- test basics
fprintf('  test ISREGULARGRID helper: basic test 1...');

[X,Y] = meshgrid(1:2,1:3);
X = X(:); Y = Y(:);
rootCoords = [X,Y];

passed = isRegularGrid(rootCoords);
if passed; fprintf(' passed!\n'); success = success + 1;
else fprintf(2,' error!\n'); failures = failures + 1;
end

clear passed X Y rootCoords;


% --- test basics
fprintf('  test ISREGULARGRID helper: basic test 2...');

[X,Y] = meshgrid(1:2,1:3);
X = X(:); Y = Y(:);
rootCoords = [X,Y];
rootCoords(end+1,:) = [0 0];

passed = ~isRegularGrid(rootCoords);
if passed; fprintf(' passed!\n'); success = success + 1;
else fprintf(2,' error!\n'); failures = failures + 1;
end

clear passed X Y rootCoords;

% --- test assertions

rootCoords = [0 0 0; 1 1 1];

fprintf('  test ISREGULARGRID helper: assertion 1...');
passed = 1;
try
  isRegularGrid(rootCoords);
  passed = 0; % should never be reached
catch ME
  if ~strcmp(ME.identifier,'isRegularGrid:dimCheck'); passed = 0; end
end
if passed; fprintf(' passed!\n'); success = success + 1;
else fprintf(2,' error!\n'); failures = failures + 1;
end

clear passed rootCoords;


% ======================================================================
% sortCoords
% ======================================================================
% --- test basics
fprintf('  test SORTCOORDS helper: basic test 1...');

coords   = [1 1; 0 0; 1 0];
expected = [0 0; 1 1; 1 0];
passed = (sum(sortCoords(coords) - expected) == 0);
if passed; fprintf(' passed!\n'); success = success + 1;
else fprintf(2,' error!\n'); failures = failures + 1;
end

clear passed coords expected;


% --- test basics
fprintf('  test SORTCOORDS helper: basic test 2...');

coords   = [1 1; 0 0; 1 0];
expected = [0 0; 1 0; 1 1];
passed = (sum(sortCoords(coords) - expected) == 0);
if passed; fprintf(' passed!\n'); success = success + 1;
else fprintf(2,' error!\n'); failures = failures + 1;
end

clear passed coords expected;

% --- test basics
fprintf('  test SORTCOORDS helper: basic test 3...');

coords   = [1 1; 0 0; 1 0];
expected = [0 0; 1 0; 1 1];
result   = sortCoords(sortCoords(coords),2);

passed = (sum(expected - result) == 0);
if passed; fprintf(' passed!\n'); success = success + 1;
else fprintf(2,' error!\n'); failures = failures + 1;
end

clear passed coords;


% --- test assertions

coords = [1 1 1; 0 0 0];

fprintf('  test SORTCOORDS helper: assertion 1...');
passed = 1;
try
  sortCoords(coords);
  passed = 0; % should never be reached
catch ME
  if ~strcmp(ME.identifier,'sortCoords:dimCheck'); passed = 0; end
end
if passed; fprintf(' passed!\n'); success = success + 1;
else fprintf(2,' error!\n'); failures = failures + 1;
end

clear passed coords;

% --- test assertions

coords = [1 1 ; 0 0];

fprintf('  test SORTCOORDS helper: assertion 2...');
passed = 1;
try
  sortCoords(coords,3);
  passed = 0; % should never be reached
catch ME
  if ~strcmp(ME.identifier,'sortCoords:column'); passed = 0; end
end
try
  sortCoords(coords,0);
  passed = 0; % should never be reached
catch ME
  if ~strcmp(ME.identifier,'sortCoords:column'); passed = 0; end
end
if passed; fprintf(' passed!\n'); success = success + 1;
else fprintf(2,' error!\n'); failures = failures + 1;
end

clear passed coords;


% ======================================================================
% isParent
% ======================================================================
domains =  [0 1 5 2 6;
            0 5 9 6 10;
            0 2 6 3 7;
            3 6 10 7 11;
            3 3 7 4 8;
            1 7 11 8 12];

% --- test basics
fprintf('  test ISPARENT helper: basic test 1...');

passed = isParent(3,domains);
if passed; fprintf(' passed!\n'); success = success + 1;
else fprintf(2,' error!\n'); failures = failures + 1;
end

clear passed;

% --- test basics
fprintf('  test ISPARENT helper: basic test 2...');

result = isParent([1 3 2],domains);
expected = [1 1 0];
passed = (sum(result - expected) == 0);
if passed; fprintf(' passed!\n'); success = success + 1;
else fprintf(2,' error!\n'); failures = failures + 1;
end

clear passed result expected;

% --- test basics
fprintf('  test ISPARENT helper: basic test 3...');

result = isParent([3 1; 2 2],domains);
expected = [1 1 ; 0 0];
passed = (sum(result - expected) == 0);
if passed; fprintf(' passed!\n'); success = success + 1;
else fprintf(2,' error!\n'); failures = failures + 1;
end

clear passed;

% --- test assertions
fprintf('  test ISPARENT helper: assertion 1...');
passed = 1;
try
  isParent(1,domains(:,1:4));
  passed = 0; % should never be reached
catch ME
  if ~strcmp(ME.identifier,'isParent:dimCheck'); passed = 0; end
end
if passed; fprintf(' passed!\n'); success = success + 1;
else fprintf(2,' error!\n'); failures = failures + 1;
end

clear passed;

% --- test assertions
fprintf('  test ISPARENT helper: assertion 2...');
passed = 1;
try
  isParent(0,domains);
  passed = 0; % should never be reached
catch ME
  if ~strcmp(ME.identifier,'isParent:i'); passed = 0; end
end
if passed; fprintf(' passed!\n'); success = success + 1;
else fprintf(2,' error!\n'); failures = failures + 1;
end

clear passed domains;


% ======================================================================
% isSemi
% ======================================================================
[X,Y] = meshgrid([3 10 7],[4 8 1 2]);
rootCoords = [X(:),Y(:)]; clear X Y;
rootCoords(end+1:end+5,:) = [8.5 2;
                           7   3;
                           8.5 3;
                           10  3;
                           8.5 4];
domains    = [0 3 11 4 12;
              0 11 7 12 8;
              0 4 12 1 9;
              0 12 8 9 5;
              0 1 9 2 10;
              0 9 5 10 6;
              4 12 13 14 15;
              4 13 8 15 16;
              4 14 15 9 17;
              4 15 16 17 5;
              2 11 7 12 13;
              2 11 7 13 8;
              3 4 12 1 14;
              3 4 14 1 9;
              6 9 17 10 6;
              6 17 5 10 6];

% --- test basics
fprintf('  test ISSEMI helper: meta test (8 tests)...');
toTest         = [1 11 12 13 14 15 16 9];
expectedReturn = [0  1  1  1  1  1  1 0];
expectedPos    = {'notSemi','SOUTH','SOUTH','WEST','WEST','NORTH','NORTH','notSemi'};
[result,pos]   = isSemi(toTest,domains,rootCoords);

passed = (result == expectedReturn).*cellfun(@(p,q) strcmp(p,q),pos,expectedPos);

if sum(passed) == length(passed); 
  fprintf(' passed!\n'); 
  success = success + sum(passed);
else 
  fprintf(2,' %i errors!\n',length(passed)-sum(passed)); 
  failures = failures + (length(passed) - sum(passed)); 
  success = success+sum(passed);
end

clear passed expectedReturn toTest expectedPos result pos domains rootCoords;


% ======================================================================
% findRoot
% ======================================================================
[X,Y] = meshgrid([3 10 7],[4 8 1 2]);
rootCoords = [X(:),Y(:)]; clear X Y;
rootCoords(end+1:end+5,:) = [8.5 2;
                           7   3;
                           8.5 3;
                           10  3;
                           8.5 4];
                           

% --- test assertions
fprintf('  test FINDROOT helper: assertion 1...');
passed = 1;
try
  findRoot(1,2,[1 2 3]);
  passed = 0; % should never be reached
catch ME
  if ~strcmp(ME.identifier,'findRoot:dimCheck'); passed = 0; end
end
if passed; fprintf(' passed!\n'); success = success + 1;
else fprintf(2,' error!\n'); failures = failures + 1;
end

clear passed domains;

% --- test assertions
fprintf('  test FINDROOT helper: assertion 2...');
passed = 1;
try
  findRoot(1,[1 2],[1 2]);
  passed = 0; % should never be reached
catch ME
  if ~strcmp(ME.identifier,'findRoot:XY'); passed = 0; end
end
if passed; fprintf(' passed!\n'); success = success + 1;
else fprintf(2,' error!\n'); failures = failures + 1;
end

clear passed domains;


% ======================================================================
% computeLimits
% ======================================================================
[X,Y] = meshgrid([3 10 7],[4 8 1 2]);
rootCoords = [X(:),Y(:)]; clear X Y;
rootCoords(end+1:end+5,:) = [8.5 2;
                           7   3;
                           8.5 3;
                           10  3;
                           8.5 4];
domains    = [0 3 11 4 12;
              0 11 7 12 8;
              0 4 12 1 9;
              0 12 8 9 5;
              0 1 9 2 10;
              0 9 5 10 6;
              4 12 13 14 15;
              4 13 8 15 16;
              4 14 15 9 17;
              4 15 16 17 5;
              2 11 7 12 13;
              2 11 7 13 8;
              3 4 12 1 14;
              3 4 14 1 9;
              6 9 17 10 6;
              6 17 5 10 6];


% --- test assertions
fprintf('  test COMPUTELIMITS helper: assertion 1...');
passed = 1;
try
  computeLimits(5,ones(4,5),ones(3,2));
  passed = 0; % should never be reached
catch ME
  if ~strcmp(ME.identifier,'computeLimits:i'); passed = 0; end
end
try
  computeLimits(0,ones(4,5),ones(3,2));
  passed = 0; % should never be reached
catch ME
  if ~strcmp(ME.identifier,'computeLimits:i'); passed = 0; end
end
if passed; fprintf(' passed!\n'); success = success + 1;
else fprintf(2,' error!\n'); failures = failures + 1;
end

clear passed domains;

% --- test assertions
fprintf('  test COMPUTELIMITS helper: assertion 2...');
passed = 1;
try
  computeLimits(5,ones(7,4),ones(3,2));
  passed = 0; % should never be reached
catch ME
  if ~strcmp(ME.identifier,'computeLimits:dimCheck'); passed = 0; end
end
if passed; fprintf(' passed!\n'); success = success + 1;
else fprintf(2,' error!\n'); failures = failures + 1;
end

clear passed domains;

% --- test assertions
fprintf('  test COMPUTELIMITS helper: assertion 3...');
passed = 1;
try
  computeLimits(5,ones(7,5),ones(3,1));
  passed = 0; % should never be reached
catch ME
  if ~strcmp(ME.identifier,'computeLimits:dimCheck'); passed = 0; end
end
if passed; fprintf(' passed!\n'); success = success + 1;
else fprintf(2,' error!\n'); failures = failures + 1;
end

clear passed domains;


% ======================================================================
% findNeighbors
% ======================================================================
[X,Y] = meshgrid([3 10 7],[4 8 1 2]);
rootCoords = [X(:),Y(:)]; clear X Y;
rootCoords(end+1:end+5,:) = [8.5 2;
                           7   3;
                           8.5 3;
                           10  3;
                           8.5 4];
domains    = [0 3 11 4 12;
              0 11 7 12 8;
              0 4 12 1 9;
              0 12 8 9 5;
              0 1 9 2 10;
              0 9 5 10 6;
              4 12 13 14 15;
              4 13 8 15 16;
              4 14 15 9 17;
              4 15 16 17 5;
              2 11 7 12 13;
              2 11 7 13 8;
              3 4 12 1 14;
              3 4 14 1 9;
              6 9 17 10 6;
              6 17 5 10 6];

% ======================================================================



nbTests = success + failures;
fprintf('\n  Number of tests:%i\n',nbTests);
fprintf('   -- %i success(es) (%1.0f %%)\n',success,100*success/nbTests);
fprintf('   -- %i failure(s) (%1.0f %%)\n\n',failures,100*failures/nbTests);

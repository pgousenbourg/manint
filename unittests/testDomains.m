% Unit tests: Helpers
% 
% Files tested:
%   * makeDomain.m
%   * findDomain.m
%   * refineDomain.m
% 
% This file is part of the project "bezierfitting" with B. Wirth from
% uni-muenster and PY. Gousenbourger from UCLouvain.
% 
% Author: Pierre-Yves Gousenbourger.
% Version: Jan. 09, 2020
% log: Dec. 17, 2019 - PYG
%        First version
%      Jan. 09, 2020 - PYG
%        Version based on testHelpers.m

close all;

addpath(genpath([pwd,'/../methods']));


disp('Unit tests on thes');
success  = 0;
failures = 0;

% Visualization
if exist('visualization') ~= 1 % (1 means that it is a variable)
  rep = input('Do you want to activate the visualization tests? [y/n]','s');
  if strcmp(rep,'y')
    visualization = 1;
    disp('');
    disp('Visualization tests activated.');
    disp('');
  elseif strcmp(rep,'n')
    visualization = 0;
    disp('');
    disp('Visualization tests deactivated.');
    disp('');
  else
    warning('I didn''t understand the response. No visualization activated.')
    disp('');
    visualization = 0;
  end
end

% ======================================================================
% makeDomain
% ======================================================================
% --- test basics
fprintf('  test MAKEDOMAIN: basic test 1...');
  
[X,Y] = meshgrid(1:3,1:4);
X = X(:); Y = Y(:); X = [X,Y];
Y = makeDomain(X);
expected = [0 1 5 2 6;
            0 5 9 6 10;
            0 2 6 3 7;
            0 6 10 7 11;
            0 3 7 4 8;
            0 7 11 8 12];

passed = (sum(Y-expected) == 0);
if passed; fprintf(' passed!\n'); success = success + 1;
else fprintf(2,' error!\n'); failures = failures + 1;
end

clear passed X Y expected;


% --- test basics
fprintf('  test MAKEDOMAIN: basic test 2...');
  
[X,Y] = meshgrid(linspace(0,1,3),linspace(1,2,4));
X = X(:); Y = Y(:); X = [X,Y];
Y = makeDomain(X);
expected = [0 1 5 2 6;
            0 5 9 6 10;
            0 2 6 3 7;
            0 6 10 7 11;
            0 3 7 4 8;
            0 7 11 8 12];

passed = (sum(Y-expected) == 0);
if passed; fprintf(' passed!\n'); success = success + 1;
else fprintf(2,' error!\n'); failures = failures + 1;
end

clear passed X Y expected;

% --- test basics
fprintf('  test MAKEDOMAIN: shuffled entries...');
  
[X,Y] = meshgrid([3 10 7],[4 8 1 2]);
X = X(:); Y = Y(:); X = [X,Y];
Y = makeDomain(X);
expected = [0 3 11 4 12;
            0 11 7 12 8;
            0 4 12 1 9;
            0 12 8 9 5;
            0 1 9 2 10;
            0 9 5 10 6];

passed = (sum(Y-expected) == 0);
if passed; fprintf(' passed!\n'); success = success + 1;
else fprintf(2,' error!\n'); failures = failures + 1;
end

clear passed X Y expected;

% --- test assertions
fprintf('  test MAKEDOMAIN: assertion 1...');
passed = 1;
try
  makeDomain([0 0 0; 1 1 1]);
  passed = 0; % should never be reached
catch ME
  if ~strcmp(ME.identifier,'surfaceFitting:dimCheck'); passed = 0; end
end
if passed; fprintf(' passed!\n'); success = success + 1;
else fprintf(2,' error!\n'); failures = failures + 1;
end

clear passed;


% --- test assertions
fprintf('  test MAKEDOMAIN: assertion 2...');
passed = 1;
try
  makeDomain([0 0; 1 1; 0 1]);
  passed = 0; % should never be reached
catch ME
  if ~strcmp(ME.identifier,'surfaceFitting:regularGrid'); passed = 0; end
end
if passed; fprintf(' passed!\n'); success = success + 1;
else fprintf(2,' error!\n'); failures = failures + 1;
end

clear passed;




% ======================================================================
% findDomain
% ======================================================================
% --- test basics
fprintf('  test FINDDOMAIN: basic test 1...');
  
[X,Y] = meshgrid(1:3,1:4);
rootCoords = [X(:),Y(:)];
domains    = [0 1 5 2 6;
              0 5 9 6 10;
              0 2 6 3 7;
              0 6 10 7 11;
              0 3 7 4 8;
              0 7 11 8 12];

passed = (findDomain(2.2,2.5,domains,rootCoords) == 4);
if passed; fprintf(' passed!\n'); success = success + 1;
else fprintf(2,' error!\n'); failures = failures + 1;
end

clear passed X Y rootCoords domains;


% --- test basics
fprintf('  test FINDDOMAIN: basic test 2...');
  
[X,Y] = meshgrid(1:3,1:4);
rootCoords = [X(:),Y(:)];
domains    = [0 1 5 2 6;
              0 5 9 6 10;
              0 2 6 3 7;
              0 6 10 7 11;
              0 3 7 4 8;
              0 7 11 8 12];

passed = (findDomain(1,1,domains,rootCoords) == 1);
if passed; fprintf(' passed!\n'); success = success + 1;
else fprintf(2,' error!\n'); failures = failures + 1;
end

clear passed X Y rootCoords domains;



% --- test matrix input
fprintf('  test FINDDOMAIN: matrix input...');
  
[X,Y] = meshgrid(1:3,1:4);
rootCoords = [X(:),Y(:)];
domains    = [0 1 5 2 6;
              0 5 9 6 10;
              0 2 6 3 7;
              0 6 10 7 11;
              0 3 7 4 8;
              0 7 11 8 12];

result   = findDomain([1 2.2; 3 1.5],[1 2.5; 3.1 3.5],domains,rootCoords);
expected = [1 4 ; 6 5];

passed = (sum(sum(result-expected)) == 0);

if passed; fprintf(' passed!\n'); success = success + 1;
else fprintf(2,' error!\n'); failures = failures + 1;
end

clear passed X Y rootCoords domains result expected;


% --- test complex
fprintf('  test FINDDOMAIN: complex test 1...');
  
[X,Y] = meshgrid([3 10 7],[4 8 1 2]);
rootCoords = [X(:),Y(:)];
domains    = [0 3 11 4 12;
              0 11 7 12 8;
              0 4 12 1 9;
              0 12 8 9 5;
              0 1 9 2 10;
              0 9 5 10 6];

passed = (findDomain(7.5,3,domains,rootCoords) == 4);

if passed; fprintf(' passed!\n'); success = success + 1;
else fprintf(2,' error!\n'); failures = failures + 1;
end

clear passed X Y rootCoords domains;


% --- test complex (refined domains)
fprintf('  test FINDDOMAIN: complex test 2 (refined)...');
  
[X,Y] = meshgrid([3 10 7],[4 8 1 2]);
rootCoords = [X(:),Y(:)];
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
              4 15 16 17 5];

passed = (findDomain(8,3.5,domains,rootCoords) == 9);
if passed; fprintf(' passed!\n'); success = success + 1;
else fprintf(2,' error!\n'); failures = failures + 1;
end

clear passed X Y rootCoords domains;

% --- test complex (refined domains with semi-domains)
fprintf('  test FINDDOMAIN: complex test 3 (refined with semis - 32 tests)...');
  
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

% 12 tests inside the domains
X(1:12) = repmat([4 8 9]',4,1);
Y(1:12) = sort(repmat([1.5,2.5,3.5,7]',3,1));
expected(1:12) = [1 11 12 13 7 8 14 9 10 5 15 16]';

% 17 tests at the corners
X(12+1:12+17) = [3 7 10 3 7 8.5 10 7 8.5 10 3 7 8.5 10 3 7 10]';
Y(12+1:12+17) = [1 1  1 2 2 2    2 3 3    3 4 4 4    4 8 8  8]';
expected(12+1:12+17) = [1 1 12 1 1 11 12 13 7 8 14 14 9 10 5 5 16]';

% 3 tests on the borders of the semis
X(end+1:end+3) = [8.5 4 8.5]';
Y(end+1:end+3) = [1.5 3 7.5]';
expected(end+1:end+3) = [11 13 15]';

X = X(:);
Y = Y(:);
expected = expected(:);

results = findDomain(X,Y,domains,rootCoords);
passed  = (results == expected);

if sum(passed) == length(X); 
  fprintf(' passed!\n'); 
  success = success + sum(passed);
else 
  fprintf(2,' %i errors!\n',length(X)-sum(passed)); 
  failures = failures + (length(X) - sum(passed)); 
  success = success+sum(passed);
  fprintf(2,'    (X,Y) values of errors:')
    errorsIdx = find(passed == 0);
    dC = [X,Y];
    XY = dC(errorsIdx,:)
end

clear passed X Y rootCoords domains results expected;


% --- test complex (refined domains with semi-domains - west excepted)
fprintf('  test FINDDOMAIN: complex test 4 (refined with semis - 32 tests)...');
  
[X,Y] = meshgrid([3 10 7],[4 8 1 2]);
rootCoords = [X(:),Y(:)]; clear X Y;
rootCoords(end+1:end+5,:) = [5 2;
                             3 3;
                             5 3;
                             7 3;
                             5 4];
domains    = [0 3 11 4 12;
              0 11 7 12 8;
              0 4 12 1 9;
              0 12 8 9 5;
              0 1 9 2 10;
              0 9 5 10 6;
              3 4 13 14 15;
              3 13 12 15 16;
              3 14 15 1 17;
              3 15 16 17 9;
              1 3 11 4 13;
              1 3 11 13 12;
              4 12 8 16 5;
              4 16 8 9 5;
              5 1 17 2 10;
              5 17 9 2 10];

% 12 tests inside the domains
X(1:12) = repmat([4 6 8]',4,1);
Y(1:12) = sort(repmat([1.5 2.5 3.5 7]',3,1));
expected(1:12) = [11 12 2 7 8 13 9 10 14 15 16 6]';

% 17 tests at the corners
X(12+1:12+17) = [3 7 10 3 5 7 10 3 5 7 3 5 7 10 3 7 10]';
Y(12+1:12+17) = [1 1  1 2 2 2  2 3 3 3 4 4 4  4 8 8  8]';
expected(12+1:12+17) = [11 12 2 11 11 12 2 7 7 8 9 9 10 14 15 16 6]';

% 3 tests on the borders of the semis
X(end+1:end+3) = [5   8 5]';
Y(end+1:end+3) = [1.5 3 7]';
expected(end+1:end+3) = [11 13 15]';

X = X(:);
Y = Y(:);
expected = expected(:);

results = findDomain(X,Y,domains,rootCoords);
passed  = (results == expected);

if sum(passed) == length(X); 
  fprintf(' passed!\n'); 
  success = success + sum(passed);
else 
  fprintf(2,' %i errors!\n',length(X)-sum(passed)); 
  failures = failures + (length(X) - sum(passed)); 
  success = success+sum(passed);
  fprintf(2,'    (X,Y) values of errors:')
    errorsIdx = find(passed == 0);
    dC = [X,Y];
    XY = dC(errorsIdx,:)
end

clear passed X Y rootCoords domains results expected;


% --- test assertions (rootCoords)
fprintf('  test FINDDOMAIN: assertion 1 (rootCoords)...');
passed = 1;
try
  findDomain(1,1,[],[0 0 0]);
  passed = 0; % should never be reached
catch ME
  if ~strcmp(ME.identifier,'findDomain:dimCheck'); passed = 0; end
end
try
  findDomain(1,1,[],0);
  passed = 0; % should never be reached
catch ME
  if ~strcmp(ME.identifier,'findDomain:dimCheck'); passed = 0; end
end
if passed; fprintf(' passed!\n'); success = success + 1;
else fprintf(2,' error!\n'); failures = failures + 1;
end

clear passed;


% --- test assertions (domains)
fprintf('  test FINDDOMAIN: assertion 2 (domains)...');
passed = 1;
try
  findDomain(1,1,0,[0 0]);
  passed = 0; % should never be reached
catch ME
  if ~strcmp(ME.identifier,'findDomain:dimCheck'); passed = 0; end
end
try
  findDomain(1,1,[0 0 0 0 0 0],[0 0]);
  passed = 0; % should never be reached
catch ME
  if ~strcmp(ME.identifier,'findDomain:dimCheck'); passed = 0; end
end
if passed; fprintf(' passed!\n'); success = success + 1;
else fprintf(2,' error!\n'); failures = failures + 1;
end

clear passed;

% --- test assertions (XY)
fprintf('  test FINDDOMAIN: assertion 3 (size XY)...');
passed = 1;
try
  findDomain(1,[1 2],ones(1,5),zeros(1,2));
  passed = 0; % should never be reached
catch ME
  if ~strcmp(ME.identifier,'findDomain:dimCheck'); passed = 0; end
end
try
  findDomain([1 2 3],[1 2],ones(1,5),zeros(1,2));
  passed = 0; % should never be reached
catch ME
  if ~strcmp(ME.identifier,'findDomain:dimCheck'); passed = 0; end
end
if passed; fprintf(' passed!\n'); success = success + 1;
else fprintf(2,' error!\n'); failures = failures + 1;
end

clear passed;


% --- test assertions (outOfRange)
fprintf('  test FINDDOMAIN: assertion 4 (outOfRange)...');
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
passed = 1;
try
  findDomain(0,2,domains,rootCoords);
  passed = 0; % should never be reached
catch ME
  if ~strcmp(ME.identifier,'findDomain:outOfRange'); passed = 0; end
end
try
  findDomain(12,2,domains,rootCoords);
  passed = 0; % should never be reached
catch ME
  if ~strcmp(ME.identifier,'findDomain:outOfRange'); passed = 0; end
end
try
  findDomain(5,10,domains,rootCoords);
  passed = 0; % should never be reached
catch ME
  if ~strcmp(ME.identifier,'findDomain:outOfRange'); passed = 0; end
end
try
  findDomain(5,0,domains,rootCoords);
  passed = 0; % should never be reached
catch ME
  if ~strcmp(ME.identifier,'findDomain:outOfRange'); passed = 0; end
end
if passed; fprintf(' passed!\n'); success = success + 1;
else fprintf(2,' error!\n'); failures = failures + 1;
end

clear passed X Y rootCoords domains;

% ======================================================================
% refineDomain
% ======================================================================
% --- test basics
fprintf('  test REFINEDOMAIN: basic test 1 (classic)...');

[X,Y] = meshgrid(0:4,0:3);
rootCoords = [X(:),Y(:)]; clear X Y;
domains = makeDomain(rootCoords);

if visualization
  figure;
  plotDomain(domains,rootCoords);
  title('REFINEDOMAIN: initialization of the domain');
end

[domains,rootCoords] = refineDomain(6,domains,rootCoords);
expectedDomains = [0     1     5     2     6;
                   0     5     9     6    10;
                   0     9    13    10    14;
                   0    13    17    14    18;
                   0     2     6     3     7;
                   0     6    10     7    11;
                   0    10    14    11    15;
                   0    14    18    15    19;
                   0     3     7     4     8;
                   0     7    11     8    12;
                   0    11    15    12    16;
                   0    15    19    16    20;
                   6     6    21    22    23;
                   6    21    10    23    24;
                   6    22    23     7    25;
                   6    23    24    25    11;
                   2     5     9     6    21;
                   2     5     9    21    10;
                   5     2     6     3    22;
                   5     2    22     3     7;
                   7    10    14    24    15;
                   7    24    14    11    15;
                   10     7    25     8    12;
                   10    25    11     8    12];

expectedRootCoords = [0    0;
                      0    1;
                      0    2;
                      0    3;
                      1    0;
                      1    1;
                      1    2;
                      1    3;
                      2    0;
                      2    1;
                      2    2;
                      2    3;
                      3    0;
                      3    1;
                      3    2;
                      3    3;
                      4    0;
                      4    1;
                      4    2;
                      4    3;
                      1.5  1;
                      1    1.5;
                      1.5  1.5;
                      2    1.5;
                      1.5  2];

if visualization
  figure;
  plotDomain(domains,rootCoords);
  title('REFINEDOMAIN: test 1 - classic');
end
tol =1e-8;
passed = (abs(sum(sum(domains - expectedDomains))) <= tol);
passed = passed*(abs(sum(sum(rootCoords - expectedRootCoords))) <= tol);

if passed; fprintf(' passed!\n'); success = success + 1;
else fprintf(2,' error!\n'); failures = failures + 1;
end

clear passed tol expectedDomains expectedRootCoords;

% Analyse the duplicates in domain
fprintf('  test REFINEDOMAIN: basic test 1 (classic) - domain duplicates...');
passed = 1;
for i = 1:size(domains,1)
  passed = passed*(length(domains(i,2:5)) - length(unique(domains(i,2:5))) == 0);
end
if passed; fprintf(' passed!\n'); success = success + 1;
else fprintf(2,' error!\n'); failures = failures + 1;
end

clear passed i;

% --- test basics 2
fprintf('  test REFINEDOMAIN: basic test 2 (neighbor of semi)...');

[domains,rootCoords] = refineDomain(1,domains,rootCoords);

expectedDomains = [0     1     5     2     6;
                   0     5     9     6    10;
                   0     9    13    10    14;
                   0    13    17    14    18;
                   0     2     6     3     7;
                   0     6    10     7    11;
                   0    10    14    11    15;
                   0    14    18    15    19;
                   0     3     7     4     8;
                   0     7    11     8    12;
                   0    11    15    12    16;
                   0    15    19    16    20;
                   6     6    21    22    23;
                   6    21    10    23    24;
                   6    22    23     7    25;
                   6    23    24    25    11;
                   2     5    26    27    28;
                   2    26     9    28    29;
                   5     2    30    31    32;
                   5    30     6    32    22;
                   7    10    14    24    15;
                   7    24    14    11    15;
                  10     7    25     8    12;
                  10    25    11     8    12;
                   2    27    28     6    21;
                   2    28    29    21    10;
                   1     1    34    35    36;
                   1    34     5    36    27;
                   3     9    13    29    14;
                   3    29    13    10    14;
                   5    31    32     3    33;
                   5    32    22    33     7;
                   1    35    36     2    30;
                   1    36    27    30     6;
                   9     3    33     4     8;
                   9    33     7     4     8];

expectedRootCoords = [0    0;
                      0    1;
                      0    2;
                      0    3;
                      1    0;
                      1    1;
                      1    2;
                      1    3;
                      2    0;
                      2    1;
                      2    2;
                      2    3;
                      3    0;
                      3    1;
                      3    2;
                      3    3;
                      4    0;
                      4    1;
                      4    2;
                      4    3;
                      1.5  1;
                      1    1.5;
                      1.5  1.5;
                      2    1.5;
                      1.5  2;
                      1.5  0;
                      1    0.5;
                      1.5  0.5;
                      2    0.5;
                      0.5  1;
                      0    1.5;
                      0.5  1.5;
                      0.5  2;
                      0.5  0;
                      0    0.5;
                      0.5  0.5];

if visualization
  figure;
  plotDomain(domains,rootCoords);
  title('REFINEDOMAIN: test 2 - neighbor of semis');
end

tol =1e-8;
passed = (abs(sum(sum(domains - expectedDomains))) <= tol);
passed = passed*(abs(sum(sum(rootCoords - expectedRootCoords))) <= tol);
if passed; fprintf(' passed!\n'); success = success + 1;
else fprintf(2,' error!\n'); failures = failures + 1;
end

clear passed tol expectedDomains expectedRootCoords;

% Analyse the duplicates in domain
fprintf('  test REFINEDOMAIN: basic test 2 (neighbor of semi) - domain duplicates...');
passed = 1;
for i = 1:size(domains,1)
  passed = passed*(length(domains(i,2:5)) - length(unique(domains(i,2:5))) == 0);
end
if passed; fprintf(' passed!\n'); success = success + 1;
else fprintf(2,' error!\n'); failures = failures + 1;
end

clear passed i;

% --- test basics 3
fprintf('  test REFINEDOMAIN: basic test 3 (semis)...');

[domains,rootCoords] = refineDomain([2.5,0.75],domains,rootCoords);

expectedDomains = [0     1     5     2     6;
                   0     5     9     6    10;
                   0     9    13    10    14;
                   0    13    17    14    18;
                   0     2     6     3     7;
                   0     6    10     7    11;
                   0    10    14    11    15;
                   0    14    18    15    19;
                   0     3     7     4     8;
                   0     7    11     8    12;
                   0    11    15    12    16;
                   0    15    19    16    20;
                   6     6    21    22    23;
                   6    21    10    23    24;
                   6    22    23     7    25;
                   6    23    24    25    11;
                   2     5    26    27    28;
                   2    26     9    28    29;
                   5     2    30    31    32;
                   5    30     6    32    22;
                   7    10    40    24    41;
                   7    40    14    41    42;
                  10     7    25     8    12;
                  10    25    11     8    12;
                   2    27    28     6    21;
                   2    28    29    21    10;
                   1     1    34    35    36;
                   1    34     5    36    27;
                   3     9    37    29    38;
                   3    37    13    38    39;
                   5    31    32     3    33;
                   5    32    22    33     7;
                   1    35    36     2    30;
                   1    36    27    30     6;
                   9     3    33     4     8;
                   9    33     7     4     8;
                   3    29    38    10    40;
                   3    38    39    40    14;
                   4    13    17    39    18;
                   4    39    17    14    18;
                   7    24    41    11    43;
                   7    41    42    43    15;
                   8    14    18    42    19;
                   8    42    18    15    19;
                  11    11    43    12    16;
                  11    43    15    12    16];

expectedRootCoords = [0    0;
                      0    1;
                      0    2;
                      0    3;
                      1    0;
                      1    1;
                      1    2;
                      1    3;
                      2    0;
                      2    1;
                      2    2;
                      2    3;
                      3    0;
                      3    1;
                      3    2;
                      3    3;
                      4    0;
                      4    1;
                      4    2;
                      4    3;
                      1.5  1;
                      1    1.5;
                      1.5  1.5;
                      2    1.5;
                      1.5  2;
                      1.5  0;
                      1    0.5;
                      1.5  0.5;
                      2    0.5;
                      0.5  1;
                      0    1.5;
                      0.5  1.5;
                      0.5  2;
                      0.5  0;
                      0    0.5;
                      0.5  0.5;
                      2.5  0;
                      2.5  0.5;
                      3    0.5;
                      2.5  1;
                      2.5  1.5;
                      3    1.5;
                      2.5  2];
                      
if visualization
  figure;
  plotDomain(domains,rootCoords);
  title('REFINEDOMAIN: test 3 - semis');
end

tol = 1e-8;
passed = (abs(sum(sum(domains - expectedDomains))) <= tol);
passed = passed*(abs(sum(sum(rootCoords - expectedRootCoords))) <= tol);
if passed; fprintf(' passed!\n'); success = success + 1;
else fprintf(2,' error!\n'); failures = failures + 1;
end

clear passed tol expectedDomains expectedRootCoords;

% Analyse the duplicates in domain
fprintf('  test REFINEDOMAIN: basic test 3 (semis) - domain duplicates...');
passed = 1;
for i = 1:size(domains,1)
  passed = passed*(length(domains(i,2:5)) - length(unique(domains(i,2:5))) == 0);
end
if passed; fprintf(' passed!\n'); success = success + 1;
else fprintf(2,' error!\n'); failures = failures + 1;
end

clear passed i;


% --- test basics 4
fprintf('  test REFINEDOMAIN: stress test 4 (cascade)...');

[domains,rootCoords] = refineDomain([2.75,1.75],domains,rootCoords);

expectedDomains = [0     1     5     2     6;
                   0     5     9     6    10;
                   0     9    13    10    14;
                   0    13    17    14    18;
                   0     2     6     3     7;
                   0     6    10     7    11;
                   0    10    14    11    15;
                   0    14    18    15    19;
                   0     3     7     4     8;
                   0     7    11     8    12;
                   0    11    15    12    16;
                   0    15    19    16    20;
                   6     6    21    22    23;
                   6    21    10    23    24;
                   6    22    23     7    25;
                   6    23    24    25    11;
                   2     5    26    27    28;
                   2    26     9    28    29;
                   5     2    30    31    32;
                   5    30     6    32    22;
                   7    10    40    24    41;
                   7    40    14    41    42;
                  10     7    25    55    56;
                  10    25    11    56    51;
                   2    27    28     6    21;
                   2    28    29    21    10;
                   1     1    34    35    36;
                   1    34     5    36    27;
                   3     9    37    29    38;
                   3    37    13    38    39;
                   5    31    32     3    33;
                   5    32    22    33     7;
                   1    35    36     2    30;
                   1    36    27    30     6;
                   9     3    33    58    59;
                   9    33     7    59    55;
                   3    29    38    10    40;
                   3    38    39    40    14;
                   4    13    48    39    49;
                   4    48    17    49    50;
                   7    24    41    11    43;
                   7    41    42    43    15;
                   8    14    44    42    45;
                   8    44    18    45    46;
                  11    11    43    51    52;
                  11    43    15    52    53;
                   8    42    45    15    47;
                   8    45    46    47    19;
                   4    39    49    14    44;
                   4    49    50    44    18;
                  12    15    47    53    61;
                  12    47    19    61    62;
                  11    51    52    12    54;
                  11    52    53    54    16;
                  10    55    56     8    57;
                  10    56    51    57    12;
                   9    58    59     4    60;
                   9    59    55    60     8;
                  12    53    61    16    63;
                  12    61    62    63    20;
                  42    41    64    65    66;
                  42    64    42    66    67;
                  42    65    66    43    68;
                  42    66    67    68    15;
                  22    40    14    41    64;
                  22    40    14    64    42;
                  41    24    41    11    65;
                  41    24    65    11    43;
                  47    42    45    67    47;
                  47    67    45    15    47;
                  46    43    68    52    53;
                  46    68    15    52    53];

expectedRootCoords = [0    0;
                      0    1;
                      0    2;
                      0    3;
                      1    0;
                      1    1;
                      1    2;
                      1    3;
                      2    0;
                      2    1;
                      2    2;
                      2    3;
                      3    0;
                      3    1;
                      3    2;
                      3    3;
                      4    0;
                      4    1;
                      4    2;
                      4    3;
                      1.5  1;
                      1    1.5;
                      1.5  1.5;
                      2    1.5;
                      1.5  2;
                      1.5  0;
                      1    0.5;
                      1.5  0.5;
                      2    0.5;
                      0.5  1;
                      0    1.5;
                      0.5  1.5;
                      0.5  2;
                      0.5  0;
                      0    0.5;
                      0.5  0.5;
                      2.5  0;
                      2.5  0.5;
                      3    0.5;
                      2.5  1;
                      2.5  1.5;
                      3    1.5;
                      2.5  2;
                      3.5  1;
                      3.5  1.5;
                      4    1.5;
                      3.5  2;
                      3.5  0
                      3.5  0.5;
                      4    0.5;
                      2    2.5;
                      2.5  2.5;
                      3    2.5;
                      2.5  3;
                      1    2.5;
                      1.5  2.5;
                      1.5  3;
                      0    2.5;
                      0.5  2.5;
                      0.5  3;
                      3.5  2.5;
                      4    2.5;
                      3.5  3;
                      2.75 1.5;
                      2.5  1.75;
                      2.75 1.75;
                      3    1.75;
                      2.75 2];

if visualization
  figure;
  plotDomain(domains,rootCoords);
  title('REFINEDOMAIN: test 4 - cascade');
end

tol = 1e-8;
passed = (abs(sum(sum(domains - expectedDomains))) <= tol);
passed = passed*(abs(sum(sum(rootCoords - expectedRootCoords))) <= tol);
if passed; fprintf(' passed!\n'); success = success + 1;
else fprintf(2,' error!\n'); failures = failures + 1;
end

clear passed tol expectedDomains expectedRootCoords;

% Analyse the duplicates in domain
fprintf('  test REFINEDOMAIN: stress test 4 (cascade) - domain duplicates...');
passed = 1;
for i = 1:size(domains,1)
  passed = passed*(length(domains(i,2:5)) - length(unique(domains(i,2:5))) == 0);
end
if passed; fprintf(' passed!\n'); success = success + 1;
else fprintf(2,' error!\n'); failures = failures + 1;
end

clear passed i domains rootCoords;

% --- test assertions (rootCoords)
fprintf('  test REFINEDOMAIN: assertion 1 (rootCoords)...');
passed = 1;
try
  refineDomain(1,[],[0 0 0]);
  passed = 0; % should never be reached
catch ME
  if ~strcmp(ME.identifier,'refineDomain:dimCheck'); passed = 0; end
end
try
  refineDomain(1,[],0);
  passed = 0; % should never be reached
catch ME
  if ~strcmp(ME.identifier,'refineDomain:dimCheck'); passed = 0; end
end
if passed; fprintf(' passed!\n'); success = success + 1;
else fprintf(2,' error!\n'); failures = failures + 1;
end

clear passed;


% --- test assertions (domains)
fprintf('  test REFINEDOMAIN: assertion 2 (domains)...');
passed = 1;
try
  refineDomain(1,0,[0 0]);
  passed = 0; % should never be reached
catch ME
  if ~strcmp(ME.identifier,'refineDomain:dimCheck'); passed = 0; end
end
try
  refineDomain(1,[0 0 0 0 0 0],[0 0]);
  passed = 0; % should never be reached
catch ME
  if ~strcmp(ME.identifier,'refineDomain:dimCheck'); passed = 0; end
end
if passed; fprintf(' passed!\n'); success = success + 1;
else fprintf(2,' error!\n'); failures = failures + 1;
end

clear passed;

% --- test assertions (I)
fprintf('  test REFINEDOMAIN: assertion 3 (size XY)...');
passed = 1;
try
  refineDomain(0,ones(1,5),zeros(1,2));
  passed = 0; % should never be reached
catch ME
  if ~strcmp(ME.identifier,'refineDomain:outOfRange'); passed = 0; end
end
try
  refineDomain(2,ones(1,5),zeros(1,2));
  passed = 0; % should never be reached
catch ME
  if ~strcmp(ME.identifier,'refineDomain:outOfRange'); passed = 0; end
end
if passed; fprintf(' passed!\n'); success = success + 1;
else fprintf(2,' error!\n'); failures = failures + 1;
end

clear passed;


% ======================================================================





nbTests = success + failures;
fprintf('\n  Number of tests:%i\n',nbTests);
fprintf('   -- %i success(es) (%1.0f %%)\n',success,100*success/nbTests);
fprintf('   -- %i failure(s) (%1.0f %%)\n\n',failures,100*failures/nbTests);

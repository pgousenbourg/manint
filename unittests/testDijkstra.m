% Unit tests: Helpers
% 
% Files tested:
%   * dijkstra.m
% 
% This file is part of the project "bezierfitting" with B. Wirth from
% uni-muenster and PY. Gousenbourger from UCLouvain.
% 
% Author: Pierre-Yves Gousenbourger.
% Version: Jan. 21, 2020
% log: Jan. 21, 2020 - PYG
%       First version

close all;

addpath(genpath([pwd,'/../methods']));
addpath(genpath([pwd,'/../manopt']));

disp('Just a silly test on Dijkstra. It'' visual...');

M = euclideanfactory(2);

domains = [0     1     5     2     6;
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

rootCoords = [0    0;
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

graph = domain2adjacency(domains,rootCoords);

% Test on dijkstra
in = 1;
out = 67;

[cost,route] = dijkstra(graph,in,out); 

% plot route
close all; 
figure; 
plotDomain(domains,rootCoords); 
hold on; 
for i = 1:length(route); 
  plot(rootCoords(route(i),1),rootCoords(route(i),2),'.r','MarkerSize',30); 
  pause(.5); 
end

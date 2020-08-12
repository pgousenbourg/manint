% Unit tests: Helpers
% 
% Files tested:
%   * rbf.m
%   * rootPos.m
%   * surfaceWeights.m
% 
% This file is part of the project "bezierfitting" with B. Wirth from
% uni-muenster and PY. Gousenbourger from UCLouvain.
% 
% Author: Pierre-Yves Gousenbourger.
% Version: Dec. 17, 2019
% log: Dec. 17, 2019 - PYG
%       First version

close all;
clear all;
clc;

successTot  = 0;
failuresTot = 0;

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


testHelpers;
  successTot = successTot + success;
  failuresTot = failuresTot + failures;

  clearvars -except successTot failuresTot visualization

testEuclideanThinPlateSpline;
  successTot = successTot + success;
  failuresTot = failuresTot + failures;

  clearvars -except successTot failuresTot visualization

testAveraging;
  successTot = successTot + success;
  failuresTot = failuresTot + failures;

  clearvars -except successTot failuresTot visualization

testBlendSurfaces;
  successTot = successTot + success;
  failuresTot = failuresTot + failures;

  clearvars -except successTot failuresTot visualization

testLiftings
  successTot = successTot + success;
  failuresTot = failuresTot + failures;

  clearvars -except successTot failuresTot visualization
  
  
testDomains
  successTot = successTot + success;
  failuresTot = failuresTot + failures;

  clearvars -except successTot failuresTot visualization

% ======================================================================
nbTestsTot = successTot + failuresTot;
fprintf('------------------------------------------------------------\n\n');
fprintf('               Total: %i tests launched.\n',nbTestsTot);
fprintf('                      %i success(es) - %1.0f %%\n',successTot,100*successTot/nbTestsTot);
fprintf('                      %i failure(s)  - %1.0f %%\n',failuresTot,100*failuresTot/nbTestsTot);

clear all

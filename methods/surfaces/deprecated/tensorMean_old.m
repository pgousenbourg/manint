function y = tensorMean(M,points,weights,varargin)
% returns the tensorized 2D-mean of points associated with weights.
%
% function Y = tensorMean(M,points)
%    returns Y, the mean of four points with equal weights.
%    The mean is computed with diadic averaging.
%
% function Y = tensorMean(M,points,weights)
%    returns Y, the mean of four points associated with given weights.
%    The mean is computed with diadic averaging.
%
% function Y = tensorMean(M,points,weights,'orientation',dim)
%    dim can be 1 or 2, and forces the diadic mean to be done
%    in the dim direction. Default: 2.
%
% inputs:  M is a manopt manifold-structure containing the exp and log
%            operators.
%          points is a 2D-structure of [p,4] manifold-valued points in
%            M. Each p sequence of points is associated with weights.
%            The 2D domain is supposed to be stacked like matlab A(:).
%          weights is a matrix (or vector) of weights associated to each
%            point.
%
% outputs: y is diadic mean of the points wrt the weights.
%
% Original author: 
% 	Pierre-Yves Gousenbourger, Dec. 18, 2019.
% Contributors: 
%
% Change log:
% 	Dec. 18, 2019 (PYG) - First version.

  % create the parser
  ip = inputParser();
  addOptional(ip,'orientation',2);

  % parse the parser as a structure.
  parse(ip, varargin{:});
  vars = ip.Results;
  
  [p,n] = size(points);
  
  % complete with missing weights
  if nargin - 2*length(varargin) == 2
    weights = ones(p,n)./(n);
  end
  
  % defense
  assert(iscell(points),'surfaceFitting:dataNotInCell','Points should be in a cell.');
  assert(size(weights,1) == p && size(weights,2) == n,'surfaceFitting:dimCheck','points and weights should have the same size')
  assert(vars.orientation == 2 || vars.orientation == 1,'tensorMean:orientation','Orientation must be 1 or 2');
  assert(n == 4,'There may only be 4 points/weights by patch.')
  
  % preprocessing
  if vars.orientation == 1
    idx = [1 3 2 4];
    points  = points(:,idx);
    weights = weights(:,idx);
  end
  
  % Dyadic operation
  y = cell(p,1);
  w = dyadicWeights(weights);
  for i = 1:p % for each set of points
    p    = dyadAv(M,points{i,1},points{i,2},weights(i,1),weights(i,2));
    q    = dyadAv(M,points{i,3},points{i,4},weights(i,3),weights(i,4));
    y{i} = dyadAv(M,p,q,w(i,1),w(i,2));
  end
end


% dyadic average
function y = dyadAv(M,x1,x2,w1,w2)
  if (w1 == 0 && w2 == 0)
    w = 0;
  else
    w = w2./(w1 + w2);
  end
  y = M.exp(x1,w.*M.log(x1,x2));
end

% weights must be given in a [px4] form.
function wNew = dyadicWeights(w);
  wNew(:,1) = (w(:,1) + w(:,2))./sum(w,2);
  wNew(:,2) = (w(:,3) + w(:,4))./sum(w,2);
end

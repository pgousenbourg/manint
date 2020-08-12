function y = blendSurfaces(M,rootPoints,tangentValues,weights,average)
% Blends four tangentValues together. 
%
% function Y = blendSurfaces(M,rootPoints,tangentValues,weights)
%    returns Y, the blended version of the four given tangentValues (in
%    the tangent space T_rootPoints(M)) and with respect to given
%    weights. By default, the averaging is done by tensorMean.
%
% function Y = blendSurfaces(M,rootPoints,tangentValues,weights,@average)
%    returns Y, the blended version of the tangentValues, with a chosen
%    averaging method, furnished as a handle function.
%
% inputs:  M, a manopt manifold-structure.
%          rootPoints, a [px4]-cell containing the rootPoints at which
%            the tangentValues are evaluated. If p = 1, then the same
%            four rootPoints will be used for every tangentValues.
%          tangentValues, a [px4] cell containing the four tangent
%            surfaces in the tangent spaces of the rootPoints.
%          weights, a [px4] matrix containing the weights used to
%            blend the tangentValues once mapped back to the manifold.
%            If p=1, then the same 4 weights will be used for every
%            tangentValues.
%          average (optional), a handleFunction that averages four
%            points with given weights. Default: tensorMean.
%            The function must be defined as
%            
%            f = @(M,points,weights) average(M,points,weights).
%
% outputs: Y, a p-cell with all the blended values in M.
%
% Original author: 
% 	Pierre-Yves Gousenbourger, Dec. 18, 2019.
% Contributors: 
%
% Change log:
% 	Dec. 18, 2019 (PYG) - First version.

  if nargin == 4
    average = @tensorMean_old;
  end
  
  % parameters
  p = size(tangentValues,1);
  
  % defense
  assert(iscell(rootPoints),'surfaceFitting:dataNotInCell','rootPoints must be in a cell');
  assert(iscell(tangentValues),'surfaceFitting:dataNotInCell','tangentValues must be in a cell');
  assert(size(rootPoints,2) == size(tangentValues,2) ...
    && size(rootPoints,2) == size(weights,2),'surfaceFitting:dimCheck',...
    'there must be as many rootPoints, tangentValues and weights to blend.');
  assert(p >= size(rootPoints,1),'surfaceFitting:dimCheck',...
    'too many lines in rootPoints compared to tangentValues');
  assert(p >= size(weights,1),'surfaceFitting:dimCheck',...
    'too many lines in weights compared to tangentValues');
  
  % generate weights and rootPoints if too few.
  if size(rootPoints,1) == 1
    rootPoints = repmat(rootPoints,p,1);
  end
  if size(weights,1) == 1
    weights = repmat(weights,p,1);
  end
  
  % Method
  %y = cell(p,1);
  points = multiExp(M,rootPoints,tangentValues);
  y      = average(M,points,weights); 
end

% maps the four tangentValues to M
%   rootPoints and tangentValues must be the same size (cell).
%   M must be a manifold.
function points = multiExp(M,rootPoints,tangentValues)
  [p,n] = size(rootPoints);
  assert(p == size(tangentValues,1) && n == size(tangentValues,2),'surfaceFitting:dimCheck','The dimensions of rootPoints and tangentValues must be the same')
  points = cell(size(rootPoints));
  for i = 1:p
  for j = 1:n
    points{i,j} = M.exp(rootPoints{i,j},tangentValues{i,j});
  end
  end
end

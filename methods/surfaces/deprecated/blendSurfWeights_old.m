function w = blendSurfWeights(t,rootCoords)
% Blending weights used in the blending technique of surfaces.
%
% function w = blendSurfWeights(t)
%    returns w, a 4-vector containing the four weights used in the
%    blending of surfaces at time t in [0,1]x[0,1]. By default, the
%    positions of the linearization points are supposed to be at the 
%    corners of this unitary patch.
%    The weights are stored in the following (graphical) manner:
%        
%         01  ---  11       ->     w(3)  ---  w(4)
%          |       |                |          |
%          |       |                |          |
%         00  ---  10              w(1)  ---  w(2)
%
% function w = blendSurfWeights(t,rootCoords)
%    returns w, a 4-vector containing the four weights used in the
%    blending of surfaces at time t in the patch defined by the 
%    rootCoords.
%    When rootCoords is not specified, it falls down naturally to
%    the unitary patch [0,1]x[0,1].
%
% EXAMPLE: 
%    >> x = 0; y = 0;
%    >> rootCoords = [0 0; 1 0; 0 1; 1 1];
%    >> blendSurfWeights(x,y,rootCoords);
%    
%      ans = 
%         
%          1  0  0  0
%    
% inputs:  t is a [px2]-vector of coordinates;
%          rootCoords is a [4x2xp] coordinates of the four corners of
%            the patch. If rootCoords is a [4x2]-patch, the same 
%            rootCoord is used for all t.
% outputs: w is a [px4]-vector of weights.
%
% Original author: 
% 	Pierre-Yves Gousenbourger, Dec. 18, 2019.
% Contributors: 
%
% Change log:
% 	Dec. 18, 2019 (PYG) - First version.

  % parameters  
  p = size(t,1);
  
  % create rootCoords if necessary
  if nargin == 1
    rootCoords = repmat([0 0; 1 0; 0 1; 1 1],1,1,size(t,1));
  end
  % rootCoords is the same for all t's.
  if (size(rootCoords,3) == 1 && p > 1)
    rootCoords = repmat(rootCoords,1,1,p);
  end

  % defense
  % on input formatting: t and rootCoords have the same size;
  assert(size(t,1) == size(rootCoords,3),'surfaceWeights:incompatibleSizes','there should be as many rootCoords as t''s');
  % on input formatting: t has max 2 entries;
  assert(size(t,2) == 2,'surfaceWeigths:incorrectT','t can have max 2 column-entries');
  
  % computation of the weights
  for i = 1:2;
  for j = 1:2;
    k  = rootPos(i,j);
    rootCoords(k,:,:);
    root = reshape(permute(rootCoords(k,:,:),[3 2 1]),[p,2]);
    x  = abs(t - root);
    
    % width of the windows
    width(:,1) = reshape(rootCoords(rootPos(2,j),1,:) - rootCoords(rootPos(1,j),1,:),p,1);
    width(:,2) = reshape(rootCoords(rootPos(i,2),2,:) - rootCoords(rootPos(i,1),2,:),p,1);
    
    %width = abs(width);
    
    assert(sum(sum((width == 0)))==0,'surfaceWeights:gridMisAlignment','The given rootCoords implied a 0-sized width.\nIs the grid turned by 90 degrees?');
    
    % weights 
    w(:,k) = prod(g(x./width),2);
  end
  end
  
end

% Auxiliary function of weights
function s = g(t)
	s = 2*t.^3 - 3*t.^2 + 1;
end

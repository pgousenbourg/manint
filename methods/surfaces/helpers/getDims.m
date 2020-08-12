function [dimGrid,dimM] = getDims(tangentValue,root)
% returns the dimGrid and dimM.
  
  sTV = size(tangentValue);
  sR = size(root);
  
  % reduce until there is no 1 at the end of the dimension
  while (sR(end) == 1 && length(sR) ~= 1)
    sR = sR(1:end-1);
  end
  dimGrid = sTV(1:end-length(sR));
  if isscalar(dimGrid)
    dimGrid = [dimGrid,1];
  elseif isempty(dimGrid)
    dimGrid = [1 1];
  end
  dimM = size(root);
end

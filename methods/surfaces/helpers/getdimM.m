function dim = getdimM(dimM)
% returns the vector of exact dimension of a vector dimM
% if dimM is a scalar, it returns a vector [dimM,1]. Otherwise, dimM is
% returned as is.
  if isscalar(dimM)
    dim = [dimM,1];
  else
    dim = dimM;
  end
end

function coeff = assembleBiStiffMat
  % Bernstein polynomials
  b = sym('[(1-s)^3,3*s*(1-s)^2,3*s^2*(1-s),s^3]');
  
  % coefficients of biharmonic stiffness matrix
  coeff = zeros(4,4,4,4);
  for i = 1:4
    for j = 1:4
      for k = 1:4
        for l = 1:4
          coeff(i,j,k,l) = int(diff(b(i),'s',2)*diff(b(k),'s',2),'s',0,1)*int(b(j)*b(l),'s',0,1)...
                          +int(b(i)*b(k),'s',0,1)*int(diff(b(j),'s',2)*diff(b(l),'s',2),'s',0,1)...
                          +2*int(diff(b(i),'s')*diff(b(k),'s'),'s',0,1)*int(diff(b(j),'s')*diff(b(l),'s'),'s',0,1);
        end
      end
    end
  end
end

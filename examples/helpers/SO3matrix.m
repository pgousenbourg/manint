% Generates SO3 matrixes.

function Q = SO3matrix(theta,phi,chi)
  if nargin == 0
    A = randn(3);
    [Q R] = qr(A);
    if (abs(det(Q) - 1) > 1)
        t = Q(1,:);
        Q(1,:) = Q(2,:);
        Q(2,:) = t;
    end
  elseif nargin == 3
    Q = rotX(theta)*rotY(phi)*rotZ(chi);
  else
    error('Wrong number of input arguments');
  end
end


function A = rotX(t)
  A =  [1 0       0     ; 
        0 cos(t) -sin(t);
        0 sin(t)  cos(t)];
end

function A = rotY(t)
  A =  [ cos(t) 0 sin(t); 
         0      1 0     ;
        -sin(t) 0 cos(t)];
end

function A = rotZ(t)
  A =  [cos(t) -sin(t) 0; 
        sin(t)  cos(t) 0;
        0       0      1];
end

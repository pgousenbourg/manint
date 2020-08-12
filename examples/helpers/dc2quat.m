function q = dc2quat(dc)
    % dc2quat quaternion direction cosine matrix angle axis
    %**********************************************************************
    %
    % dc2quat calculates the quaternion corresponding to a direction
    % cosine matrix. I believe this is based on the algorithm
    % given by A. R. Klumpp, "Singularity-Free Extraction of a
    % Quaternion from a Direction-Cosine Matrix," Journal of
    % Spacecraft and Rockets, Vol. 13, Dec. 1976, pp. 754-755.
    % Assumes input dc is orthonormal.
    %
    % Input: dc = 3x3 direction cosine matrix
    %
    % Output: q = quaternion, q(1) = scalar, q(2:4) = vector
    % Rotation sense = Successive rotations are right multiplies.
    %
    % Programmer: James Tursa
    %
    %**********************************************************************

    q = [0 0 0 0];
    tr = dc(1,1) + dc(2,2) + dc(3,3);
    ii = 0;
    nn = 0;
    q(1) = tr;
    for kk=1:3
        if( dc(kk,kk) > q(1) )
            ii = kk;
            q(1) = dc(kk,kk);
        end
    end

    tr = sqrt(1 + 2*q(1) - tr);

    order = [2 3 1];
    for mm=1:3
        kk = order(mm);
        nn = nn + 1;
        jj = 6 - nn - kk;
        x = ii * (nn - ii);
        if( x == 0)
            q(1) = (dc(jj,kk) - dc(kk,jj)) / tr;
            q(nn+1) = q(1);
        else
            q(jj+kk-ii+1) = (dc(jj,kk) + dc(kk,jj)) / tr;
        end
    end

    if( ii == 0 )
        q(1) = tr;
    else
        q(ii+1) = tr;
    end
    q(2:4) = -q(2:4);
    if( q(1) == 0 )
        q = 0.5 * q;
    else
        q = 0.5 * sign(q(1)) * q;
    end

    %\
    % MAKES QUATERNION A POSITIVE ROTATION
    %/
%    if( q(1) <= 0 )
%        q = -q;
%    end

    %\
    % NORMALIZE QUATERNION (KEEPS ROUTINE STABLE)
    %/
    q = q / norm(q);

    return

end


%------------------------------------------------------------------------

function dc = quat2dc(q)
% quat2dc quaternion direction cosine matrix angle axis
%*******************************************************************
%
% quat2dc calculates the dirction cosine matrix corresponding to a
% quaternion. Assumes input quaternion is normalized.
%
% Input: q = quaternion, q(1) = scalar, q(2:4) = vector
% Rotation sense = Successive rotations are right multiplies.
%
% Output: dc = 3x3 direction cosine matrix
%
% Programmer: James Tursa
%
%*******************************************************************

q11 = q(1)^2;
q12 = q(1)*q(2);
q13 = q(1)*q(3);
q14 = q(1)*q(4);
q22 = q(2)^2;
q23 = q(2)*q(3);
q24 = q(2)*q(4);
q33 = q(3)^2;
q34 = q(3)*q(4);
q44 = q(4)^2;

dc = zeros(3);
dc(1,1) = q11 + q22 - q33 - q44;
dc(2,1) = 2 * (q23 - q14);
dc(3,1) = 2 * (q24 + q13);
dc(1,2) = 2 * (q23 + q14);
dc(2,2) = q11 - q22 + q33 - q44;
dc(3,2) = 2 * (q34 - q12);
dc(1,3) = 2 * (q24 - q13);
dc(2,3) = 2 * (q34 + q12);
dc(3,3) = q11 - q22 - q33 + q44;

return
end

%------------------------------------------------------------------------

function [a v] = angleaxis(x)
% angleaxis quaternion direction cosine matrix angle axis
%*******************************************************************
%
% angleaxis calculates the rotation angle and rotation axis of the
% input quaternion or direction cosine matrix.
%
% Input: x = quaternion, x(1) = scalar, x(2:4) = vector
% Rotation sense = Successive rotations are right multiplies.
% Assumes x is normalized.
%
% or
%
% x = direction cosine matrix.
% Assumes x is orthonormalized.
%
% Output: a = rotation angle (radians)
% v = rotation axis (1x3 unit vector)
%
% Programmer: James Tursa
%
%*******************************************************************

if( numel(x) == 9 )
    q = dc2quat(x);
else
    q = x;
end
a = 2 * acos(q(1));
if( nargout == 2 )
    if( a == 0 )
        v = [1 0 0];
    else
        v = q(2:4)/sin(a/2);
    end
end

return
end

% Generates SO3 matrixes.

function Q = SO3matrix()
    A = randn(3);
    [Q R] = qr(A);
    if (abs(det(Q) - 1) > 1)
        t = Q(1,:);
        Q(1,:) = Q(2,:);
        Q(2,:) = t;
    end
end
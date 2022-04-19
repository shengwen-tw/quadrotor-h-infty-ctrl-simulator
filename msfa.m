%matrix sign function algorithm
function [T1, T2, fail]=msfa(H)
n = max(size(H));
I = eye(n);
%rcond(H)
%disp(min(eig(H)));
H_last = H;
H_now = [];

%solve sign(H)
while 1
    %[U, S, V] = svd(H_last);
    %inv_H_last = V / S * U';

    inv_H_last = inv(H_last);
    %disp(norm(H_last * inv_H_last))
    
    H_now = (H_last + inv_H_last) / 2;
    residual = norm(H_now - H_last);
    H_last = H_now;
    if residual < 1e-6
        break;
    end
end

Q = null(I - H_now); %solve (I - sign(H)) * X = 0
%[Q, R] = qr(H_now - I);

T1 = Q(1:n/2, 1:n/2);
T2 = Q(n/2+1:n, 1:n/2);

fail = 0;
end

%check "calculating the inertia of a real symmetric (or tridiagonal) matrix"
%on "math.stackexchange.com"
function flag=is_psd_matrix(T1, T2)
    flag = 1;
    %[P, Omega] = tridiag(T1.'*T2);
    %[L, D] = ldl_tri(Omega);
    %disp(Omega);
    %disp(D);
    [V, D] = eig(T1.'*T2);
    
    for i = 1:max(size(Omega))
        if(D(i, i) < 0)
            flag = 0;
            return;
        end
    end
end
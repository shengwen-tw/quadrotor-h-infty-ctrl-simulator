%symmetric matrix tridiagonalization
function [P, T]=tridiag(A)
    mat_size = max(size(A));
    T = A;
    P = eye(mat_size, mat_size);
    for i = 1:(mat_size - 2)
        sub_mat_size = mat_size - i;
        
        e = zeros(sub_mat_size, 1);
        e(1) = 1;
        x = T(i + 1:end, i);
        norm_x = norm(x, 2);
        a = sign(e.'*x) * norm_x;
        ut_x = sqrt(1/2*norm_x*(norm_x + abs(x(1))));
        u = (x + a*e) / (2*ut_x);
        H = eye(sub_mat_size, sub_mat_size) - 2*((u*u.')/(u.'*u));
        P_i = eye(mat_size, mat_size);
        P_i(i+1:end, i+1:end) = H;
        
        T(i, i+1:end) = T(i, i+1:end) * H;
        T(i+1:end, i) = H * T(i+1:end, i);
        T(i+1:end,i+1:end) = H * T(i+1:end,i+1:end) * H.';
        
        P = P_i * P;
    end
end
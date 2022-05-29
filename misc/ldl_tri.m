%LDL.' decomposition for tridiagonal matrix
function [L, D]=ldl_tri(T)
    mat_size = max(size(T));
    D = T;
    L = eye(mat_size, mat_size);
    for i = 1:(mat_size - 1)
        L_i = eye(mat_size, mat_size);
        L_i(i + 1, i) = -D(i+1, i) / D(i, i);
        D = L_i * D * L_i.';
        L = L_i * L;
    end
end
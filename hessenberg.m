%matrix hessenberg form transformation
function [Q, H]=hessenberg(A)
mat_size = max(size(A));
H = A;
Q = eye(mat_size, mat_size);
%select column
for i = 1:(mat_size - 2)    
    %eliminate row
    for k = 0:(mat_size - i - 2)
        %Givens rotation on column
        ii = mat_size - k - 1;
        jj = mat_size - k;
        
        c =  H(ii, i) / sqrt(H(ii, i)^2 + H(jj, i)^2);
        s = -H(jj, i) / sqrt(H(ii, i)^2 + H(jj, i)^2);
        
        G = eye(mat_size, mat_size);
        G(ii, ii) = c;
        G(jj, jj) = c;
        G(ii, jj) = -s;
        G(jj, ii) = s;
        
        H = G * H * G.';
        
        Q = G * Q;
    end
end
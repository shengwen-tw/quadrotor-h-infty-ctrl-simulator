function [Q, S, eig_val]=qr_iteration(A, k)
n = max(size(A));
H = A;

%convert input matrix to hessenberg form
[Q0, H] = hessenberg(A);
Q = Q0.';

%QR algorithm
for i = 1:k
    [Q_i, R] = qr(H);
    H = R * Q_i;
    Q = Q * Q_i;
end
S = H;

%calculate eigenvalues
i = 1;
eig_val = [];
while i <= n
    %check if reach the end of the row or not
    if ((i+1) <= n)
        %1x1 block
        if abs(S(i + 1, i)) < 1e-6
            %1x1 block
            eig_val = [eig_val; S(i, i)]; %real eigenvalue
        else
            %2x2 block
            eig_val = [eig_val; solve2x2_eigen(S(i:i+1,i:i+1))]; %complex eigenvalue
            i = i + 1;
        end
    elseif i == n
        %last element on the diagonal
        eig_val = [eig_val; S(i, i)];
    end
    
    i = i + 1;
end
%disp(eig_val);
end

function[eig]=solve2x2_eigen(A)
eig = [0; 0];
a = A(1, 1);
b = A(1, 2);
c = A(2, 1);
d = A(2, 2);
a_plus_d = a + d;
ad = a*d;
bc = b*c;
eig(1) = ((a_plus_d)+sqrt((a_plus_d)*(a_plus_d)-4*(ad-bc)))/2;
eig(2) = ((a_plus_d)-sqrt((a_plus_d)*(a_plus_d)-4*(ad-bc)))/2;
end

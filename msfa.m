function [T1, T2, fail]=msfa(H)
n = max(size(H));
I = eye(n);
%H = H - 0.8*sqrt(-1)*I;
%disp(min(eig(H)));
H_last = balance(H);
H_now = [];
while 1
    H_now = (H_last + inv(H_last)) / 2;
    residual = norm(H_now - H_last);
    H_last = H_now;
    if(residual < 1e-10)
        break;
    end
end
[Q, R] = qr(H_now - I);
T1 = Q(1:n/2, 1:n/2);
T2 = Q(n/2+1:n, 1:n/2);
fail = 0;
end
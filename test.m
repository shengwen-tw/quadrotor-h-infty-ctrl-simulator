%quadrotor parameters
g = 9.8;
m = 1;
Ix = 0.01466;
Iy = 0.01466;
Iz = 0.02848;

%system dynamics
s_phi = sin(0);
c_phi = cos(0);
s_theta = sin(0);
c_theta = cos(0);
s_psi = sin(0);
c_psi = cos(0);
t_theta = tan(0);
sec_theta = sec(0);
p = 0;
q = 0;
r = 0;
u = 0;
v = 0;
w = 0;

a1 = [-r*s_phi*t_theta + q*c_phi*t_theta ...
    r*(c_phi*sec_theta^2 + q*s_phi*sec_theta^2) ...
    0 1 s_phi*t_theta c_phi*t_theta 0 0 0 0 0 0];
a2 = [(-q*s_phi - r*c_phi) 0 0 0 c_phi -s_phi 0 0 0 0 0 0];
a3 = [-r*s_phi/c_theta + q*c_phi/c_theta ...
    r*c_phi*sec_theta*t_theta + q*s_phi*sec_theta*t_theta ...
    0 0 s_phi/c_theta c_phi/c_theta 0 0 0 0 0 0];
a4 = [0 0 0 0 (Iy-Iz)/Ix*r (Iy-Iz)/Ix*q 0 0 0 0 0 0];
a5 = [0 0 0 (Iz-Ix)/Iy*r 0 (Iz-Ix)/Iy*p 0 0 0 0 0 0];
a6 = [0 0 0 (Ix-Iy)/Iz*q (Ix-Iy)/Iz*p 0 0 0 0 0 0 0];
a7 = [0 -g*c_theta 0 0 -w v 0 r -q 0 0 0];
a8 = [g*c_phi*c_theta -g*s_phi*s_theta 0 w 0 -u -r 0 p 0 0 0];
a9 = [-g*c_theta*s_phi -g*s_theta*c_phi 0 -v u 0 q -p 0 0 0 0];
a10 = [w*(c_phi*s_psi - s_phi*c_psi*s_theta) + v*(s_phi*s_psi + c_psi*c_phi*s_theta) ...
    w*(c_phi*c_psi*c_theta) + v*(c_psi*s_phi*c_theta) - u*(c_psi*s_theta) ...
    w*(s_phi*c_psi - c_phi*s_psi*s_theta) - v*(c_phi*c_psi - c_phi*c_psi*s_theta) + u*(c_theta*c_psi) ...
    0 0 0 c_psi*c_theta (-c_phi*s_psi + c_psi*s_phi*s_theta) (s_phi*s_psi + c_phi*c_psi*s_theta) 0 0 0];
a11 = [v*(-s_phi*c_psi + c_phi*s_psi*s_theta) - w*(c_psi*c_phi + s_phi*s_psi*s_theta) ...
    v*(s_phi*s_psi*c_theta) + w*(c_phi*s_psi*c_theta) - u*(s_theta*s_psi) ...
    v*(-c_phi*s_psi + s_phi*c_psi*s_theta) + w*(s_psi*s_phi + c_phi*c_psi*s_theta) + u*(c_theta*c_psi) ...
    0 0 0 c_theta*s_psi (c_phi*c_psi + s_phi*s_psi*s_theta) (-c_psi*s_phi + c_phi*s_psi*s_theta) 0 0 0];
a12 = [-w*s_phi*c_theta + v*c_theta*c_phi ...
    -w*c_phi*s_theta - u*c_theta - v*s_theta*s_phi ...
    0 0 0 0 -s_theta c_theta*s_phi c_phi*c_theta 0 0 0];
A = [a1; a2; a3; a4; a5; a6; a7; a8; a9; a10; a11; a12];

C1 = zeros(16, 12);
C1(1, 1) = 200;    %roll
C1(2, 2) = 200;    %pitch
C1(3, 3) = 500;    %yaw
C1(4, 4) = 10;     %roll rate
C1(5, 5) = 10;     %pitch rate
C1(6, 6) = 20;     %yaw rate
C1(7, 7) = 700;    %vx
C1(8, 8) = 700;    %vy
C1(9, 9) = 500;    %vz
C1(10, 10) = 2000; %x
C1(11, 11) = 2000; %y
C1(12, 12) = 1000; %z

B1 = [0    0   0   0   0   0;
      0    0   0   0   0   0;
      0    0   0   0   0   0;
      0    0   0  1/Ix 0   0;
      0    0   0   0  1/Iy 0;
      0    0   0   0   0  1/Iz;
    -1/m   0   0   0   0   0;
      0 -1/m   0   0   0   0;
      0    0 -1/m  0   0   0;
      0    0   0   0   0   0;
      0    0   0   0   0   0;
      0    0   0   0   0   0;];
  
B2 = [0   0   0   0;
      0   0   0   0;
      0   0   0   0;
      0  1/Ix 0   0;
      0   0  1/Iy 0;
      0   0   0  1/Iz;
      0   0   0   0;
      0   0   0   0;
    -1/m  0   0   0;
      0   0   0   0;
      0   0   0   0;
      0   0   0   0];
  
D12 = [0 0 0 0;
       0 0 0 0;
       0 0 0 0;
       0 0 0 0;
       0 0 0 0;
       0 0 0 0;
       0 0 0 0;
       0 0 0 0;
       0 0 0 0;
       0 0 0 0;
       0 0 0 0;
       0 0 0 0;
       1 0 0 0;
       0 1 0 0;
       0 0 1 0;
       0 0 0 1];

%H-infinity synthesis
if 0
[gamma, X] = hinf_syn(A, B1, B2, C1, 0);
disp(sprintf("gamma: %d\ngamma(db):%d", gamma, mag2db(gamma)))
%
C0_hat = -B2.' * X;
sys=ss(A + B2*C0_hat, B1, C1 + D12*C0_hat, 0);
sigma(sys, ss(gamma));
end

if 0
A = [-1/2 -1/2; 0 3/10];
C1 = [1 0; 0 1; 0 0; 0 0];
B1 = C1.';
D12 = [0 0; 0 0; 1 0; 0 1];
D21 = D12.';
B2= [1/5 3/2; -3/2 0];
C2 = [-7/5 -7/10; -1/2 -7/10];
[gamma, X] = hinf_syn(A, B1, B2, C1, 0);
%
C0_hat = -B2.' * X;
sys=ss(A + B2*C0_hat, B1, C1 + D12*C0_hat, 0);
sigma(sys, ss(gamma));
end

if 0
A = [-1];
B1 = [1 0];
B2 = [1];
C1 = [1; 0];
C2 = [1];
D12 = [0; 1];
D21 = [1 0]
[gamma, X] = hinf_syn(A, B1, B2, C1, 0);
disp('r_infimum')
disp(gamma)
%
C0_hat = -B2.' * X;
sys=ss(A + B2*C0_hat, B1, C1 + D12*C0_hat, 0);
sigma(sys, ss(gamma));
end

%test symmetric matrix tridiagonalization
if 0
A = [ 4 1 -2  2;
      1 2  0  1;
     -2 0  3 -2;
      2 1 -2 -1];
[P, T] = tridiag(A);
answer = [ 4  -3    0      0;
          -3 10/3 -5/3     0;
           0 -5/3 -33/25 68/75;
           0   0   68/75 149/75];
disp(answer);
disp(T);
disp(P * A * P.') %T = P*A*Pt
disp(P);
end

%test diagonalization for tridiagonalized matrix
if 0
T = [ 4  -3    0      0;
     -3 10/3 -5/3     0;
      0 -5/3 -33/25 68/75;
      0   0   68/75 149/75]
  
[L, D] = ldl_tri(T);
disp(D);
disp(L*T*L.')
%disp(L);
end

%test hessenberg form transformation
if 0
A =[16     2     3    13
     5    11    10     8
     9     7     6    12
     4    14    15     1];
%A = magic(5);
[Q, H] = hessenberg(A);
disp(A);
disp(H);
disp(Q*A*Q.')
end

%test qr iteration
if 0
A = [2.00000    3.00000    1.00000    0.50000    4.00000
     4.00000    5.00000    7.00000    0.10000    1.00000
     5.00000    3.00000    6.00000   19.20000    9.00000
     1.00000    4.00000    1.00000    4.00000    7.00000
     3.00000    1.00000    6.00000    2.00000    6.00000];
disp("MATLAB:");
S = schur(A)
eig(S)
disp("Mine:");
[Q, S, eig_val] = qr_iteration(A, 10000);
disp(S);
disp(eig_val);
disp(Q.' * A * Q)
end
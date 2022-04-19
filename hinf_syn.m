function [gamma, X, X_norm]=hinf_syn(A, B1, B2, C1, D)
eps = 1e-6;
care_norm_eps = 1e-4;
iteration = 0;

%calculate lower bound gamma
At = A.';
B2t = B2.';
H = B2 * B2t;
G = C1.' * C1;
[Z, fail] = care_sda(At, H, G); %a fixed-point iteration algorithm for solving CARE
%z_norm = norm(A*Z + Z*At - Z*G*Z + H);
%disp(z_norm)
B1t = B1.';
C1t = C1.';
gamma_l = hinf_norm((A - Z*C1t*C1).', C1t, B1t, 0);
%disp('gamma_lb:');
%disp(gamma_l);

%approximate an upper bound gamma:
%this is not the exact upper bound of the gamma. user should check
%the value is feasible before using it.
rho_max = eigs(B1.'*B1, 1,'lm'); %largest singular value of the B1
rho_0 = eigs(B2.'*B2, 1,'sm');   %smallest singular value of the B2
gamma_u = rho_max / rho_0;
%disp('gamma_ub');
%disp(gamma_u);

optimal_norm = 1e10;
optimal_gamma = gamma_u;
optimal_X = [];

%bisection method for optimal gamma and CARE solution
while 1
    iteration = iteration + 1;
    
    gamma = (gamma_l + gamma_u) / 2;
    
    %construct Hamiltonian matrix
    H = [   A,     1/(gamma*gamma).*(B1*B1t) - B2*B2t;
         -C1t*C1,                   -A.'            ];
    
    if(has_pure_img_eigen(H) == 0)
        %[T1, T2, fail] = msfa(H);
        %X = T2 / T1;
        
        H = C1.'*C1;
        G = B2*B2.' - (1 / (gamma*gamma) .* B1*B1.');
        [X, fail] = care_sda(A, H, G);
                
        %check if the solution is stable
        if(fail == 0)
            X_norm = norm(At*X + X*A - X*G*X + H);
            if(X_norm < care_norm_eps)
                %found a good gamma and CARE solution
                optimal_norm = X_norm;
                optimal_gamma = gamma;
                optimal_X = X;
                %
                gamma_u = gamma;
            else
                gamma_l = gamma;
            end
        else
            gamma_l = gamma;
        end
    else
        gamma_l = gamma;
    end
    
    %bisection value is converged
    if (gamma_u - gamma_l) < eps
        if(optimal_norm < care_norm_eps)
            %H = [   A,     1/(gamma*gamma).*(B1*B1t) - B2*B2t;
            %     -C1t*C1,                   -A.'            ];
            %[T1, T2, fail] = msfa(H);
            %X = T2 / T1;

            gamma = optimal_gamma;
            X = optimal_X;
            return;
        else
            %no solution found, return H2 control solution
            H = C1.'*C1;
            G = B2*B2.';
            X = care_sda(A, 0, H, G);
            X_norm = norm(At*X + X*A - X*G*X + H);
            %gamma = 0; %just convenient for plotting
            return;
        end
    end
end
end

function [T1, T2]=get_stable_invariant_subspace(H)
[V, D] = eig(H);
[m, n] = size(H);

V_stable = [];
D_stable = zeros(m/2, n/2);
j = 1;
for i = 1:m
    if(real(D(i, i)) < -1e-6)
        V_stable = [V_stable, V(:, i)];
        D_stable(j, j) = D(i, i);
        j = j + 1;
    end
end

T1 = V_stable(1:m/2, :);
T2 = V_stable(m/2+1:m, :);

%disp(diag(D_stable))
%disp(norm(H*V_new - V_new*D_new));
end

function retval=has_pure_img_eigen(H)
retval = 0;
[V, D] = eig(H);

for i= 1 : max(size(H))
    if(~isreal(D(i, i)) && abs(real(D(i, i))) < 1e-2)
        retval = 1;
        return;
    end
end
end

%check "calculating the inertia of a real symmetric (or tridiagonal) matrix"
%on "math.stackexchange.com"
function flag=is_psd_matrix(T1, T2)
    flag = 1;
    [P, Omega] = tridiag(T1.'*T2);
    [L, D] = ldl_tri(Omega);
    
    %disp(Omega);
    %disp(D);
    
    %check inertia of the matrix
    for i = 1:max(size(Omega))
        %if exisits negative element of the D matrix,
        %X is not P.S.D
        if(D(i, i) < 0)
            flag = 0;
            return;
        end
    end
end
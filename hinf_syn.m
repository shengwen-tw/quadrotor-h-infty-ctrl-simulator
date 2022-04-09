function [gamma, X]=hinf_syn(A, B1, B2, C1, D)
    eps = 1e-6;
    iteration = 0;

    %calculate lower bound gamma
    At = A.';
    B2t = B2.';
    H = B2 * B2t;
    G = C1.' * C1;
    Z = care_sda(At, 0, H, G); %a fixed-point iteration algorithm for solving CARE
    %z_norm = norm(A*Z + Z*At - Z*G*Z + H);
    %disp(z_norm)
    B1t = B1.';
    C1t = C1.';
    %gamma_l = hinf_norm((A - Z*C1t*C1).', C1t, B1t, 0);
    gamma_l = hinf_norm((A - Z*C1t*C1).', C1t, B1t, 0) + 100;
    
    %approximate an upper bound gamma:
    %this is not the exact upper bound of the gamma. user should check
    %the value is feasible before using it.
    rho_max = eigs(B1.'*B1, 1,'lm'); %largest singular value of the B1
    rho_0 = eigs(B2.'*B2, 1,'sm');   %smallest singular value of the B2
    gamma_u = rho_max / rho_0;
    %gamma_u = 1e10;
    
    %bisection and secant method start
    while(abs(gamma_u - gamma_l) > eps)
        iteration = iteration + 1;
        
        %bisection searching
        gamma = (gamma_l + gamma_u) / 2;
    
        %construct Hamiltonian matrix
        H = [   A,         1/(gamma*gamma).*(B1*B1t) - B2*B2t;
             -C1t*C1,                       -A.'            ];
    
        if(has_pure_img_eigen(H) == 0)
            [T1, T2] = get_stable_invariant_subspace(H);      
            X = real(T2 * inv(T1)); %we extract the real part of X only since the
                                    %imaginary part may remained by numerical errors
            
            if(is_psd_matrix(X) == 1)
               gamma_u = gamma; %decrease gamma upper bound
            else
               gamma_l = gamma; %increase gamma lower bound
            end         
        else
            gamma_l = gamma; %increase gamma lower bound
        end    
    end
    
    %disp(iteration);
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
    
    %disp(norm(H*V_new - V_new*D_new));
end

function retval=has_pure_img_eigen(H)
    retval = 0;
    [V, D] = eig(H);
        
    for i= 1 : max(size(H))
        if(~isreal(D(i, i)) && abs(real(D(i, i))) < 1e-6)
            retval = 1;
            return;
        end
    end
end

function retval=is_psd_matrix(A)
    retval = 1;
    [V, D] = eig(A);
    for i= 1 : max(size(A))
        if(real(D(i, i)) < 0)
            retval = 0;
            return;
        end
    end
end
function [gamma, X, X_norm]=hinf_syn(A, B1, B2, C1, D)
    eps = 1e-6;
    care_norm_eps = 1e-4;
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
            H = C1.'*C1;
            G = B2*B2.' - (1 / (gamma*gamma) .* B1*B1.');
            [X, sda_converged] = modified_sda(A, 0, H, G);
                        
            %check if the solution is stable
            if(sda_converged == 1)
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
                gamma = optimal_gamma;
                X = optimal_X;
                return;
            else                
                %no solution found, return H2 control solution
                H = C1.'*C1;
                G = B2*B2.';
                X = care_sda(A, 0, H, G);
                X_norm = norm(At*X + X*A - X*G*X + H);
                gamma = 0; %just convenient for plotting
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
        if(~isreal(D(i, i)) && abs(real(D(i, i))) < 1e-3)
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

function [X, state] = modified_sda(A, B, H, G)
state = 1;
state_dimension = length(A);

r = 2.4; %SDA's author suggested the value between 2.1~2.6
I = eye(state_dimension);
A_r = A - (r*I);

iteration_times = 0;

%solve CARE with SDA
A_hat_last = I + 2*r*inv(A_r + G*transpose(inv(A_r))*H);
G_hat_last = 2*r*inv(A_r)*G*inv(A_r.' + H*inv(A_r)*G);
H_hat_last = 2*r*inv(A_r.' + H*inv(A_r)*G)*H*inv(A_r);
norm_H_now = norm(H_hat_last);

while 1
    iteration_times = iteration_times + 1;
    
    %reduce redundent calculation by pre-calculating repeated terms
    inv_I_plus_H_G = inv(I + (H_hat_last * G_hat_last));
    A_hat_last_t = A_hat_last.';
    
    %update
    A_hat_new = A_hat_last * inv(I + G_hat_last * H_hat_last) * A_hat_last;
    G_hat_new = G_hat_last + (A_hat_last * G_hat_last * inv_I_plus_H_G * A_hat_last_t);
    H_Hat_new = H_hat_last + (A_hat_last_t * inv_I_plus_H_G * H_hat_last * A_hat_last);
        
    %matrix norms
    norm_H_last = norm_H_now;
    norm_H_now = norm(H_Hat_new);
    
    %prepare next iteration
    A_hat_last = A_hat_new;
    G_hat_last = G_hat_new;
    H_hat_last = H_Hat_new;
    
    %diverge
    if (norm_H_now - norm_H_last) < 0
        X = [];
        state = 0;
        return;
    end
    
    %stop iteration if converged
    if abs(norm_H_now - norm_H_last) < 1e-16
        break;
    end
    
    %disp(iteration_times);
end

X = H_Hat_new;
end
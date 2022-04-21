%a bisection method to sythesis the H-infinity control solution with
%minimal gamma using the structure-preserving doubling algorithm
function [gamma, X, ric_residual]=hinf_syn(A, B1, B2, C1, D)     
eps = 1e-6;
residual_eps = 1e-7;

%record of total iteration numbers
iteration = 0;

%common variables
At = A.';
B1B1t = B1*B1.';
B2B2t = B2*B2.';
C1tC1 = C1.'*C1;

%calculate lower bound gamma
H = B2B2t;
G = C1tC1;
[Z, fail] = care_sda(At, H, G);
gamma_l = hinf_norm((A - Z*C1tC1).', C1.', B1.', 0);
%disp(norm(A*Z + Z*At - Z*G*Z + H));
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

optimal_residual = 1e10;
optimal_gamma = gamma_u;
optimal_X = [];

%bisection searching for optimal gamma and CARE solution
while 1
    iteration = iteration + 1;
    
    gamma = (gamma_l + gamma_u) / 2;
    
    H = C1tC1;
    G = B2B2t - (B1B1t / (gamma*gamma));
    
    %construct Hamiltonian matrix
    Ham = [ A, -G;
           -H, -At];
    Ham = balance(Ham);
    
    if(has_pure_img_eigen(Ham) == 0)
        [X, fail] = care_sda(A, H, G);
        
        %if SDA manage to produce a solution
        if(fail == 0)
            %check if the solution is stable
            ric_residual = norm(At*X + X*A - X*G*X + H);
            if(ric_residual < residual_eps)
                %found a good gamma and CARE solution
                optimal_residual = ric_residual;
                optimal_gamma = gamma;
                optimal_X = X;
                %
                gamma_u = gamma; %decrease the upper bound
            else
                gamma_l = gamma; %increase the lower bound
            end
        else
            gamma_l = gamma; %increase the lower bound
        end
    else
        gamma_l = gamma; %increase the lower bound
    end
    
    %bisection searching is converged
    if (gamma_u - gamma_l) < eps
        if(optimal_residual < residual_eps)
            %accept the solution
            gamma = optimal_gamma;
            X = optimal_X;
            %disp(iteration);
            return;
        else
            %no solution found, return H2 control solution
            G = B2B2t;
            X = care_sda(A, H, G);
            ric_residual = norm(At*X + X*A - X*G*X + H);
            gamma = 0; %just for convenient
            return;
        end
    end
end
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

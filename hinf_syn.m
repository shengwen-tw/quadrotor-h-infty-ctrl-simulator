function gamma=hinf_syn(A, B, C, D)
    eps = 1e-6;
    gamma_l = 0;
    gamma_u = 1e10;

    while(abs(gamma_u - gamma_l) > eps)
        %bisection searching
        gamma = (gamma_l + gamma_u) / 2;
    
        %construct Hamiltonian matrix
        H = [      A,         (1/gamma)*(B*B');
             -(1/gamma)*C.'*C,       -A.'    ];
    
        if(has_pure_img_eigen(H) == 0)
            get_stable_invariant_subspace(H);
            gamma_u = gamma; %decrease gamma upper bound          
        else
            gamma_l = gamma; %increase gamma lower bound
        end    
    end
end

function [T1, T2]=get_stable_invariant_subspace(H)
    [V, D] = eig(H);
    [m, n] = size(H);
    
    V_stable = [];
    D_stable = zeros(12, 12);
    j = 1;
    for i = 1:m
        if(real(D(i, i)) > 0)
            continue;
        end

        V_stable = [V_stable, V(:, i)];
        D_stable(j, j) = D(i, i);
        j = j + 1;
    end
    
    T1 = V_stable(1:m/2);
    T2 = V_stable(m/2+1: m);
    
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
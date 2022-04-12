function gamma=hinf_norm(A, B, C, D)
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
        gamma_u = gamma; %decrease gamma upper bound
    else
        gamma_l = gamma; %increase gamma lower bound
    end
end
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
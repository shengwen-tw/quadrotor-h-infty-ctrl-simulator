function [X, fail] = care_sda(A, H, G)
fail = 0;
state_dimension = length(A);

r = 2.4; %SDA's author suggested the value between 2.1~2.6
I = eye(state_dimension);
A_r = A - (r*I);

iteration_times = 0;

%solve CARE with SDA
A_hat_last = I + 2*r*inv(A_r + G*transpose(inv(A_r))*H);
G_hat_last = 2*r*(A_r\G/(A_r.' + H/A_r*G));
H_hat_last = 2*r*inv(A_r.' + H/A_r*G)*H/A_r;
norm_H_now = norm(H_hat_last);

while 1
    iteration_times = iteration_times + 1;
    
    %reduce redundent calculation by pre-calculating repeated terms
    I_plus_H_G = I + (H_hat_last * G_hat_last);
    A_hat_last_t = A_hat_last.';
    
    %update
    A_hat_new = A_hat_last * ((I + G_hat_last * H_hat_last) \ A_hat_last);
    G_hat_new = G_hat_last + (A_hat_last * G_hat_last * (I_plus_H_G \ A_hat_last_t));
    H_Hat_new = H_hat_last + (A_hat_last_t * (I_plus_H_G \ H_hat_last * A_hat_last));
    
    %matrix norms
    norm_H_last = norm_H_now;
    norm_H_now = norm(H_Hat_new);
    
    %prepare next iteration
    A_hat_last = A_hat_new;
    G_hat_last = G_hat_new;
    H_hat_last = H_Hat_new;
    
    %diverge
    %if (norm_H_now - norm_H_last) < 0
    %    X = [];
    %    fail = 1;
    %    return;
    %end
    
    %stop iteration if converged
    if abs(norm_H_now - norm_H_last) < 1e-16
        break;
    end
    
    %disp(iteration_times);
end

X = H_Hat_new;
end
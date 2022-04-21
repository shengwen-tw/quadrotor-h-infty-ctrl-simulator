function [X, fail] = care_sda(A, H, G)
fail = 0;
max_iteration = 50;
state_dimension = length(A);

r = 2.4; %SDA's author suggested the value between 2.1~2.6
I = eye(state_dimension);
A_r = A - (r*I);

iteration = 0;

%solve CARE with SDA
A_hat_last = I + 2*r*inv(A_r + G*transpose(inv(A_r))*H);
G_hat_last = 2*r*(A_r\G/(A_r.' + H/A_r*G));
H_hat_last = 2*r*inv(A_r.' + H/A_r*G)*H/A_r;
residual_now = norm(H_hat_last);

while 1
    iteration = iteration + 1;
    
    %reduce redundent calculation by pre-calculating repeated terms
    I_plus_H_G = I + (H_hat_last * G_hat_last);
    A_hat_last_t = A_hat_last.';
    
    %update
    A_hat_new = A_hat_last * ((I + G_hat_last * H_hat_last) \ A_hat_last);
    G_hat_new = G_hat_last + (A_hat_last * G_hat_last * (I_plus_H_G \ A_hat_last_t));
    H_hat_new = H_hat_last + (A_hat_last_t * (I_plus_H_G \ H_hat_last * A_hat_last));
    
    %matrix norms
    residual_last = residual_now;
    residual_now = norm(H_hat_new - H_hat_last);
    
    %prepare next iteration
    A_hat_last = A_hat_new;
    G_hat_last = G_hat_new;
    H_hat_last = H_hat_new;
    
    %diverge
    if iteration > max_iteration
        X = [];
        fail = 1;
        return;
    end
    
    %stop iteration if converged
    if residual_now < 1e-10
        break;
    end
end
%disp(iteration_times);
X = H_hat_new;
end
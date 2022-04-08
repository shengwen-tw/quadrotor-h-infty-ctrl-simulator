function quadrotor_sim
math = se3_math;

%====================%
% Simulation options %
%====================%
ITERATION_TIMES = 10000;
COMPARE_MATLAB_SDA = 1;

%================%
% UAV parameters %
%================%
uav_dynamics = dynamics;        %create uav dynamics object
uav_dynamics.dt = 0.001;        %set iteration period [sec]
uav_dynamics.mass = 1;          %set uav mass [kg]
uav_dynamics.a = [0; 0; 0];     %acceleration of uav [m/s^2], effected by applied force
uav_dynamics.v = [0; 0.5; 0];     %initial velocity of uav [m/s]
uav_dynamics.x = [0.5; 0; 0];     %initial position of uav [m]
uav_dynamics.W = [0; 0; 0];     %initial angular velocity of uav
uav_dynamics.W_dot = [0; 0; 0]; %angular acceleration of uav, effected by applied moment
uav_dynamics.f = [0; 0; 0];     %force generated by controller
uav_dynamics.M = [0; 0; 0];     %moment generated by controller
uav_dynamics.J = [0.01466 0 0;  %inertia matrix
                  0 0.01466 0;
                  0 0 0.02848];
              
%set initial attitude with Euler angles
init_attitude(1) = deg2rad(0);  %roll
init_attitude(2) = deg2rad(0);  %pitch
init_attitude(3) = deg2rad(0);  %yaw
uav_dynamics.R = math.euler_to_dcm(init_attitude(1), init_attitude(2), init_attitude(3));

quad_sim_greeting(uav_dynamics, ITERATION_TIMES, init_attitude);

%=========================================%
% Parameters of the H-infinity controller %
%=========================================%
m = uav_dynamics.mass;
g = 9.8;
Ix = uav_dynamics.J(1, 1);
Iy = uav_dynamics.J(2, 2);
Iz = uav_dynamics.J(3, 3);

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

C1 = zeros(16, 12);
C1(1, 1) = 300;    %roll
C1(2, 2) = 300;    %pitch
C1(3, 3) = 750;    %yaw
C1(4, 4) = 10;     %roll rate
C1(5, 5) = 10;     %pitch rate
C1(6, 6) = 10;     %yaw rate
C1(7, 7) = 900;    %vx
C1(8, 8) = 900;    %vy
C1(9, 9) = 900;    %vz
C1(10, 10) = 2500; %x
C1(11, 11) = 2500; %y
C1(12, 12) = 3000; %z

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
       1 0 0 0;  %f
       0 1 0 0;  %tau_x
       0 0 1 0;  %tau_y
       0 0 0 1]; %tau_z
   
C2 = eye(12, 12);
D21 = 0;

I_6x6 = eye(6, 6);
B1B1t = B1 * B1.';
B2B2t = B2 * B2.';
C1tC1 = C1.' * C1;
B2t = B2.';

%controller setpoints
xd = zeros(3, ITERATION_TIMES);
vd = zeros(3, ITERATION_TIMES);
yaw_d = zeros(1, ITERATION_TIMES);

%plot datas
time_arr = zeros(1, ITERATION_TIMES);
speed_inc_arr = zeros(1, ITERATION_TIMES);
sda_time_arr = zeros(1, ITERATION_TIMES);

matlab_x_norm_arr = zeros(1, ITERATION_TIMES);
sda_x_norm_arr = zeros(1, ITERATION_TIMES);

matlab_time_arr = zeros(1, ITERATION_TIMES);
vel_arr = zeros(3, ITERATION_TIMES);
R_arr = zeros(3, 3, ITERATION_TIMES);
euler_arr = zeros(3, ITERATION_TIMES);
pos_arr = zeros(3, ITERATION_TIMES);
W_arr = zeros(3, ITERATION_TIMES);
M_arr = zeros(3, ITERATION_TIMES);
d_arr = zeros(6, ITERATION_TIMES);

%%%%%%%%%%%%%%%%%%%%%
%   path planning   %
%%%%%%%%%%%%%%%%%%%%%
% cirular trajectory
radius = 0.5;         %[m]
circum_rate = 0.25;   %[hz], times of finished a circular trajectory per second
climb_rate = -0.05;
yaw_rate = 0.05;      %[hz], times of full rotation around z axis per second
for i = 1: ITERATION_TIMES
    %plan heading
    if i == 1
        yaw_d(1) = 0;
    else
        yaw_d(i) = yaw_d(i - 1) + (yaw_rate * uav_dynamics.dt * 2 * pi);
    end
    if yaw_d(i) > pi   %bound yaw angle between +-180 degree
        yaw_d(i) = yaw_d(i) - (2 * pi);
    end
    
    %plan position
    xd(1, i) = radius * cos(circum_rate * uav_dynamics.dt * i * pi);
    xd(2, i) = radius * sin(circum_rate * uav_dynamics.dt * i * pi);
    xd(3, i) = i * uav_dynamics.dt * climb_rate;
    
    %end height
    if(xd(3, i) <= -1)
        xd(3, i) = -1;
    end
    
    %plan velocity
    vd(1, i) = radius * -sin(circum_rate * uav_dynamics.dt * i * pi);
    vd(2, i) = radius * cos(circum_rate * uav_dynamics.dt * i * pi);
    vd(3, i) = climb_rate;
end

progress_tok = waitbar(0, 'Starting');
for i = 1: ITERATION_TIMES
    %disp(i);
    prompt = sprintf('Progress: %d %%\n(%d/%d)', floor(i/ITERATION_TIMES*100), i, ITERATION_TIMES);
    waitbar(i/ITERATION_TIMES, progress_tok, prompt);
    
    %========================%
    % Update System Dynamics %
    %========================%
    uav_dynamics = update(uav_dynamics);
    
    %=======================%
    % Quadrotor LQR Control %
    %=======================%
    
    eulers = math.dcm_to_euler(uav_dynamics.R); %get euler angles from R matrix
    v_b = uav_dynamics.R * uav_dynamics.v;      %get body frame velocity
    
    p = uav_dynamics.W(1);
    q = uav_dynamics.W(2);
    r = uav_dynamics.W(3);
    u = v_b(1);
    v = v_b(2);
    w = v_b(3);
    
    %construct A matrix
    s_phi = sin(eulers(1));
    c_phi = cos(eulers(1));
    s_theta = sin(eulers(2));
    c_theta = cos(eulers(2));
    s_psi = sin(eulers(3));
    c_psi = cos(eulers(3));
    t_theta = tan(eulers(2));
    sec_theta = sec(eulers(2));
    
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
    At = transpose(A);
    
    %====================%
    % H-infinity control %
    %====================%
    
    %construct state vector
    x = [eulers(1);
         eulers(2);
         eulers(3);
         uav_dynamics.W(1);
         uav_dynamics.W(2);
         uav_dynamics.W(3);
         v_b(1);
         v_b(2);
         v_b(3);
         uav_dynamics.x(1);
         uav_dynamics.x(2);
         uav_dynamics.x(3)];
    
    %construct desired setpoint vector
    x0 = [deg2rad(0); %desired roll
          deg2rad(0); %desired pitch
          deg2rad(0); %yaw_d(i); %desired yaw
          0;          %desired roll
          0;          %desired pitch
          0;          %desired yaw
          vd(1, i);   %desired x velocity
          vd(2, i);   %desired y velocity
          vd(3, i);   %desired z velocity
          xd(1, i);   %desired x position
          xd(2, i);   %desired y position
          xd(3, i)];  %desired z position
     
    %=========================================================%
    % solve CARE (Continuous-time Algebraic Riccati Equation) %
    %=========================================================%
    
    %H-infinity control synthesis
    %tstart = tic();
    [gamma, X_] = hinf_syn(A, B1, B2, C1, 0);
    %sda_time = toc(tstart);
    %inv_r2 = 1 / (gamma*gamma);
    %r2_B1B1t_B2B2t = -((inv_r2 .* B1B1t) - B2B2t);
    %sda_x_norm = norm(At*X_ + X_*A - X_*r2_B1B1t_B2B2t*X_ + C1tC1)
    
    %method1: SDA (Structure-Preserving Doubling Algorithm)
    inv_r2 = 1 / (gamma*gamma);
    r2_B1B1t_B2B2t = -((inv_r2 .* B1B1t) - B2B2t);
    %
    tstart = tic();
    X = care_sda(A, B2, C1tC1, r2_B1B1t_B2B2t);
    sda_time = toc(tstart);
    sda_x_norm = norm(At*X + X*A - X*r2_B1B1t_B2B2t*X + C1tC1);
        
    %method2: MATLAB
    if COMPARE_MATLAB_SDA ~= 0
        B = [B1, B2];
        Bt = B.';
        m1 = size(B1, 2);
        m2 = size(B2, 2);
        R = [-gamma^2*eye(m1) zeros(m1, m2);
                 zeros(m2,m1)      eye(m2)];
        BRBt = B * R * Bt;
        %
        tstart = tic();
        [X_, L, G_dummy] = care(A, B, C1tC1, R);
        matlab_time = toc(tstart);
        matlab_x_norm = norm(At*X_ + X_*A - X_*BRBt*X_ + C1tC1);
    
        %efficiency comparison of CARE solvers
        speed_inc = matlab_time / sda_time;
    end
    
    %calculate feedback control
    C0_hat = -B2t * X;
    u_fb = C0_hat * [x - x0];
        
    %calculate feedforward control
    gravity_ff = dot(uav_dynamics.mass .* g .* [0; 0; 1], uav_dynamics.R * [0; 0; 1]);
    u_ff = [gravity_ff; 0; 0; 0];    
    
    %obtain complete control input
    u = u_ff + u_fb;
    
    lqr_f = uav_dynamics.R * [0; 0; u(1)];
    lqr_M = [u(2); u(3); u(4)];
    
    %feed control to uav dynamics
    uav_dynamics.f = lqr_f;
    uav_dynamics.M = lqr_M;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Record datas for plotting %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    time_arr(i) = i * uav_dynamics.dt;
    if COMPARE_MATLAB_SDA ~= 0
        sda_time_arr(i) = sda_time;
        matlab_time_arr(i) = matlab_time;
        matlab_x_norm_arr(i) = matlab_x_norm;
        sda_x_norm_arr(i) = sda_x_norm;
        speed_inc_arr(i) = speed_inc;
        sda_time_arr(i) = sda_time;
    end
    vel_arr(:, i) = uav_dynamics.v;
    pos_arr(:, i) = uav_dynamics.x;
    R_arr(:, :, i) = uav_dynamics.R;
    euler_arr(:, i) = rad2deg(math.dcm_to_euler(uav_dynamics.R));
    W_arr(:, i) = rad2deg(uav_dynamics.W);
    M_arr(:, i) = uav_dynamics.M;
    d_arr(:, i) = uav_dynamics.d;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Animate the simulation result %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rigidbody_visualize([5; 5; 5], pos_arr, R_arr, ITERATION_TIMES, uav_dynamics.dt, 30);

%%%%%%%%%%%%%%%%%%%%%%%%%%
%          Plot          %
%%%%%%%%%%%%%%%%%%%%%%%%%%

%attitude (euler angles)
figure('Name', 'attitude (euler angles)');
subplot (3, 1, 1);
plot(time_arr, euler_arr(1, :));
xlabel('time [s]');
ylabel('roll [deg]');
subplot (3, 1, 2);
plot(time_arr, euler_arr(2, :));
xlabel('time [s]');
ylabel('pitch [deg]');
subplot (3, 1, 3);
plot(time_arr, euler_arr(3, :));
xlabel('time [s]');
ylabel('yaw [deg]');

%angular velocity
figure('Name', 'Angular velocity');
subplot (3, 1, 1);
plot(time_arr, W_arr(1, :));
xlabel('time [s]');
ylabel('x [deg/s]');
subplot (3, 1, 2);
plot(time_arr, W_arr(2, :));
xlabel('time [s]');
ylabel('y [deg/s]');
subplot (3, 1, 3);
plot(time_arr, W_arr(3, :));
xlabel('time [s]');
ylabel('z [deg/s]');

%velocity
figure('Name', 'velocity (NED frame)');
subplot (3, 1, 1);
plot(time_arr, vel_arr(1, :), time_arr, vd(1, :));
xlabel('time [s]');
ylabel('x [m/s]');
subplot (3, 1, 2);
plot(time_arr, vel_arr(2, :), time_arr, vd(2, :));
xlabel('time [s]');
ylabel('y [m/s]');
subplot (3, 1, 3);
plot(time_arr, -vel_arr(3, :), time_arr, -vd(3, :));
xlabel('time [s]');
ylabel('-z [m/s]');

%position
figure('Name', 'position (NED frame)');
subplot (3, 1, 1);
plot(time_arr, pos_arr(1, :), time_arr, xd(1, :));
xlabel('time [s]');
ylabel('x [m]');
subplot (3, 1, 2);
plot(time_arr, pos_arr(2, :), time_arr, xd(2, :));
xlabel('time [s]');
ylabel('y [m]');
subplot (3, 1, 3);
plot(time_arr, -pos_arr(3, :), time_arr, -xd(3, :));
xlabel('time [s]');
ylabel('-z [m]');

if COMPARE_MATLAB_SDA ~= 0
    %time cost of the CARE solvers
    figure('Name', 'time cost of CARE solvers');
    title('time cost');
    plot(time_arr, sda_time_arr, time_arr, matlab_time_arr);
    xlabel('time [s]');
    ylabel('cost [s]');
    legend('SDA CARE', 'MATLAB CARE');

    %speed improvement of SDA compare to the MATLAB CARE
    figure('Name', 'SDA / MATLAB CARE');
    plot(time_arr, speed_inc_arr);

    %precision of CARE solvers
    disp(mean(speed_inc_arr));
    figure('Name', 'precision of CARE solvers');
    title('precision (norm of CARE)');
    plot(time_arr, sda_x_norm_arr, time_arr, matlab_x_norm_arr);
    xlabel('time [s]');
    ylabel('CARE norm');
    legend('SDA CARE', 'MATLAB CARE');
end

%disturbance
figure('Name', 'disturbances');
subplot (3, 2, 1);
plot(time_arr, d_arr(1, :));
xlabel('time [s]');
ylabel('f_{wx}');
subplot (3, 2, 3);
plot(time_arr, d_arr(2, :));
xlabel('time [s]');
ylabel('f_{wy}');
subplot (3, 2, 5);
plot(time_arr, d_arr(3, :));
xlabel('time [s]');
ylabel('f_{wz}');
subplot (3, 2, 2);
plot(time_arr, d_arr(4, :));
xlabel('time [s]');
ylabel('\tau_{wx}');
subplot (3, 2, 4);
plot(time_arr, d_arr(5, :));
xlabel('time [s]');
ylabel('\tau_{wy}');
subplot (3, 2, 6);
plot(time_arr, d_arr(6, :));
xlabel('time [s]');
ylabel('-\tau_{wz}');

disp("Press any key to leave");
pause;
close all;
delete(progress_tok);
end

function quad_sim_greeting(dynamics, iteration_times, init_attitude)
roll = rad2deg(init_attitude(1));
pitch = rad2deg(init_attitude(2));
yaw = rad2deg(init_attitude(3));
disp(sprintf('Quadrotor simulation (%d iterations, dt = %dsec)', iteration_times, dynamics.dt));
disp(sprintf('Initial position: (%f, %f, %f)', dynamics.x(1), dynamics.x(2), dynamics.x(3)));
disp(sprintf('Initial attitude (euler angle): (%f, %f, %f)', roll, pitch, yaw));
disp('Start simulation...');
end

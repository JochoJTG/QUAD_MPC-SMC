clear all
clc;

%%

%Drone parameters

m = 1.121; %kg
g = 9.81; %m/s^2

Jxx = 0.0466;        % kg*m^2
Jyy = 0.0589;        % kg*m^2
Jzz = 0.0802;        % kg*m^2

J = diag([Jxx Jyy Jzz]);

%% %Z MPC parameters

N = 10;
ts_z = 0.1;

Q_z = 10;      %500
Qt_z = 0.005;     %-0.005

A_z = [
0,1;
0,0;
];
B_z = [0;1/m];
G = [0;-g];
C = [1 0];
[A_zd, B_zd] = c2d(A_z,B_z,ts_z);

[n,nu] = size(B_zd);
ny = size(C,1);

H = 0;
F1 = 0;
F2 = 0;
F3 = 0;
F4 = 0;


Phi = [];
Psi = [];
G_mat = [];

u_des = m*g;

PHI_I = A_zd;
Aux_PSI_I = B_zd;
Aux_THETA_I = eye(n,n);
auxG_i = G;
for i = 1:N
    
    PSI_I = [Aux_PSI_I, zeros(n,(N-i)*nu)];
    %THETA_I = [Aux_THETA_I, zeros(n,(N-i)*n)];
    G_i = [auxG_i,zeros(n,(N-i)*nu)];
    SumG_1 = sum(G_i(1,:)); 
    SumG_2 = sum(G_i(2,:));
    SumG = [SumG_1;SumG_2];

    Phi = [Phi PHI_I];
    Psi = [Psi;PSI_I];

    pi_n = PI_i(i,ny,N);
    pi_nu = PI_i(i,nu,N);


   
    H = H + 2*(PSI_I'*C'*Q_z*C*PSI_I + pi_nu'*Qt_z*pi_nu);
    F1 = F1 + 2*(PSI_I'*C'*Q_z*C*PHI_I); % falta multiplicar por w(k)
    F2 = F2 - 2*(PSI_I'*C'*Q_z*pi_n); % falta multiplicar por ytilde
    F3 = F3 - 2*(pi_nu'*Qt_z); % falta multiplicar por vdeseada
    %F4 = F4 + 2*(PSI_I'*C'*Q_z*C*THETA_I);
    F4 = F4 - 2*(PSI_I'*C'*Qt_z*C*SumG);
    

    PHI_I = PHI_I*A_zd;
    Aux_PSI_I = [A_zd*Aux_PSI_I B_zd];
    %Aux_THETA_I = [A_zd*Aux_THETA_I eye(n,n)];
    auxG_i = [A_zd*auxG_i G];
  
end


%% XY MPC parameters

N_xy = 15;
ts_xy = 0.1;

Q_xy = diag([10 10]);
R_xy = diag([15 15]);


A_xy = [
0,1,0,0;
0,0,0,0;
0,0,0,1;
0,0,0,0;
];
B_xy = [
0,0;
g,0;
0,0;
0,-g;
];

C_xy = [
1 0 0 0;
0 0 1 0;
];

[A_xy_d, B_xy_d] = c2d(A_xy,B_xy,ts_xy);
[n_xy, nu_xy] = size(B_xy);
ny_xy = size(C_xy,1);

H_xy = 0;
F1_xy = 0;
F2_xy = 0;
F3_xy = 0;

Phi_xy = [];
Psi_xy = [];

u_des_xy = [0;0];

PHI_xy_I = A_xy_d;
Aux_PSI_xy_I = B_xy_d;

for i = 1:N_xy
    
    PSI_xy_I = [Aux_PSI_xy_I, zeros(n_xy,(N_xy-i)*nu_xy)];
    Phi_xy = [Phi_xy PHI_xy_I];
    Psi_xy = [Psi_xy;PSI_xy_I];

    pi_n = PI_i(i,ny_xy,N_xy);
    pi_nu = PI_i(i,nu_xy,N_xy);

    F1_xy = F1_xy + 2*(PSI_xy_I'*C_xy'*Q_xy*C_xy*PHI_xy_I); % falta multiplicar por x(k)
    F2_xy = F2_xy - 2*(PSI_xy_I'*C_xy'*Q_xy*pi_n); % falta multiplicar por ytilde
    F3_xy = F3_xy - 2*(pi_nu'*R_xy); % falta multiplicar por udeseada

    if i < N
        pi_2nu = PI_i(i+1,nu_xy,N_xy);
        H_xy = H_xy + 2*(PSI_xy_I'*C_xy'*Q_xy*C_xy*PSI_xy_I + pi_nu'*R_xy*pi_nu + (pi_2nu - pi_nu)'*R_xy*(pi_2nu - pi_nu));

    else
        H_xy = H_xy + 2*(PSI_xy_I'*C_xy'*Q_xy*C_xy*PSI_xy_I + pi_nu'*R_xy*pi_nu);
    end

    PHI_xy_I = PHI_xy_I*A_xy_d;
    Aux_PSI_xy_I = [A_xy_d*Aux_PSI_xy_I B_xy_d];
  
end





%% Setup Linear Movement Controller Parameters

N_l = 15;        %Prediction horizon for the linear displacement 10
ts_l = 0.1;    

% Qy_l = [50 0; 0 50];     %1000 500
% Qu_l = [1 0; 0 1]; %75 50

Qy_l = [20 0; 0 20];     %10 10
Qu_l = [105 0; 0 105]; %1 1

ud_l = [0;0];
          
A_lc = [0 1 0 0; 0 0 0 0; 0 0 0 1; 0 0 0 0];
B_lc = [0 0; g 0; 0 0; 0 -g];
C_l = [1 0 0 0;0 0 1 0];
D_l = [0 0; 0 0];
[A_l,B_l] = c2d(A_lc,B_lc,ts_l); 

[n_l,nu_l] = size(B_l);
nr_l = size(C_l,1);

ul_max = [pi/12;pi/12];    %Setup the upper and lower input constraints (RP commands)
ul_min = [-pi/12;-pi/12];

yl_max = [2;2];
yl_min = [-2;-2];
% yl_max = [pi/8;pi/8];
% yl_min = [-pi/8;-pi/8];

delmax_l = [0.2;0.2];
delmin_l = [-0.2;-0.2];

[max_ul,min_ul] = InConstraints(ul_max,ul_min,N_l);
[H_l,F1_l,F2_l,F3_l,F4_l,phi_l,psi_l] = Cost_Funct(A_l,B_l,C_l,G,Qy_l,Qu_l,N_l,2);
psi_last = psi_l(N_l*n_l-(n_l-1):N_l*n_l,:);
[Aineq_l,G1_l,G2_l,G3_l] = Ineq_Calc(C_l,phi_l,psi_l ,N_l,nr_l,nu_l,n_l,yl_max,yl_min,delmax_l,delmin_l,B_l);

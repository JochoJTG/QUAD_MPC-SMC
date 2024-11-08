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
ts_z = 0.2;

Q_z = 500000;      %500
Qt_z = 0.0000001;     %-0.005

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

u_des = 0;

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
    F4 = F4 + 2*(PSI_I'*C'*Qt_z*C*SumG);

    PHI_I = PHI_I*A_zd;
    Aux_PSI_I = [A_zd*Aux_PSI_I B_zd];
    %Aux_THETA_I = [A_zd*Aux_THETA_I eye(n,n)];
    auxG_i = [A_zd*auxG_i G];
  
end


%% XY MPC parameters

N_xy = 10;
ts_xy = 0.2;

Q_xy = diag([10 10]);
R_xy = diag([10 10]);


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


    H_xy = H_xy + 2*(PSI_xy_I'*C_xy'*Q_xy*C_xy*PSI_xy_I + pi_nu'*R_xy*pi_nu);
    F1_xy = F1_xy + 2*(PSI_xy_I'*C_xy'*Q_xy*C_xy*PHI_xy_I); % falta multiplicar por x(k)
    F2_xy = F2_xy - 2*(PSI_xy_I'*C_xy'*Q_xy*pi_n); % falta multiplicar por ytilde
    F3_xy = F3_xy - 2*(pi_nu'*R_xy); % falta multiplicar por udeseada

    PHI_xy_I = PHI_xy_I*A_xy_d;
    Aux_PSI_xy_I = [A_xy_d*Aux_PSI_xy_I B_xy_d];
  
end





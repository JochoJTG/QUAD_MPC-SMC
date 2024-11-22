function [Aineq,G1,G2,G3] = Ineq_Calc(Cc,phi,psi,N,nr,nu,n,y_max,y_min,delmax,delmin,B)

Aineq_1 = zeros(2*N*nr,N*nu);
Aineq_2 = zeros(2*N*nu,N*nu);
G1_1 = zeros(2*N*nr,n);

G2_2 = zeros(2*N*nu,nu);
G2_2(1:nu,1:nu) = eye(nu);
G2_2(N*nu+1:N*nu+nu,1:nu) = -eye(nu);

G3_1 = zeros(2*N*nr,1);
G3_2 = zeros(2*N*nr,1);

for i = 1:N
    Aineq_1((i*nu)-(nu-1):(i*nu),1:N*nu) = Cc*PI_i(i,n,N)*psi;              %Aineq_1 Matrix
    Aineq_1((i*nu+N*nu)-(nu-1):(i*nu+N*nu),1:N*nu) = -Cc*PI_i(i,n,N)*psi;

%      Aineq_2((i*nu)-(nu-1):i*nu,(i*nu)-(nu-1):i*nu) = -Cc*psi(:,i*nu-(nu-1):i*nu);
%      Aineq_2((i*nu+N*nu)-(nu-1):(i*nu+N*nu),(i*nu)-(nu-1):i*nu) = Cc*psi(:,i*nu-(nu-1):i*nu);
%     Aineq_2((i*nu)-(nu-1):i*nu,(i*nu)-(nu-1):i*nu) = eye(nu);
%     Aineq_2((i*nu+N*nu)-(nu-1):(i*nu+N*nu),(i*nu)-(nu-1):i*nu) = -eye(nu);
    
%     Aineq_2((i*nu)-(nu-1):i*nu,(i*nu)-(nu-1):i*nu) = eye(nu);               %Aineq_2 Matrix
%     Aineq_2((i*nu+N*nu)-(nu-1):(i*nu+N*nu),(i*nu)-(nu-1):i*nu) = -eye(nu);
%     if i > 1 
%         Aineq_2((i*nu)-(nu-1):i*nu,(i*nu)-(2*nu-1):(i*nu)-(nu)) = -eye(nu);
%         Aineq_2((i*nu+N*nu)-(nu-1):(i*nu+N*nu),(i*nu)-(2*nu-1):(i*nu)-(nu)) = eye(nu);
%     end
   
    G1_1((i*nu)-(nu-1):(i*nu),1:n) = -Cc*PI_i(i,n,N)*phi;                   %G1_1 Matrix
    G1_1((i*nu+N*nu)-(nu-1):(i*nu+N*nu),1:n) = Cc*PI_i(i,n,N)*phi;

%     G1_1((i*nu)-(nu-1):(i*nu),1:n) = -Cc*phi(n*N-(n-1):n*N,:);                   %G1_1 Matrix
%     G1_1((i*nu+N*nu)-(nu-1):(i*nu+N*nu),1:n) = Cc*phi(n*N-(n-1):n*N,:);

    if nr == 2
        G3_1(i*nr,1) = y_max(nr,1);
        G3_1(i*nr-1,1) = y_max(nr-1,1);
        G3_1(N*nr+i*nr,1) = -y_min(nr,1);
        G3_1(N*nr+i*nr-1,1) = -y_min(nr-1,1);
   
%         G3_2(i*nr,1) = delmax(nr,1);
%         G3_2(i*nr-1,1) = delmax(nr-1,1);
%         G3_2(N*nr+i*nr,1) = -delmin(nr,1);
%         G3_2(N*nr+i*nr-1,1) = -delmin(nr-1,1);
        
    elseif nr == 3
        G3_1(i*nr,1) = y_max(nr,1);
        G3_1(i*nr-1,1) = y_max(nr-1,1);
        G3_1(i*nr-2,1) = y_max(nr-2,1);
        G3_1(N*nr+i*nr,1) = -y_min(nr,1);
        G3_1(N*nr+i*nr-1,1) = -y_min(nr-1,1);
        G3_1(N*nr+i*nr-2,1) = -y_min(nr-2,1);
        
%         G3_2(i*nr,1) = delmax(nr,1);
%         G3_2(i*nr-1,1) = delmax(nr-1,1);
%         G3_2(i*nr-2,1) = delmax(nr-2,1);
%         G3_2(N*nr+i*nr,1) = -delmin(nr,1);
%         G3_2(N*nr+i*nr-1,1) = -delmin(nr-1,1);
%         G3_2(N*nr+i*nr-2,1) = -delmin(nr-2,1);
        
     end
    
end

% Aineq = [Aineq_1;Aineq_2];                                                  % Complete Aineq
% G1 = [G1_1;zeros(2*N*nu,n)];                                                % Complete G1
% G2 = [zeros(2*N*nr,nu);G2_2];                                               % Complete G2
% G3 = [G3_1;G3_2];

Aineq = Aineq_1;                                                  % Complete Aineq
G1 = G1_1;                                                % Complete G1
G2 = zeros(2*N*nr,nu);                                               % Complete G2
G3 = G3_1;



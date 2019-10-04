% Robyn Woollands 2017
% Texas A&M University - Department of Aerospace Engineering
% File name     : lecture3_example9_bvp_conv.m
% Description   : Eigenvalue convergence analysis for Lecture 3 of the 
%                 JUNKINS & WOOLLANDS, Picard-Chebyshev Lecture Series.
% Date Written  : March 31, 2017
% Date Modified : March 31, 2017
%================================================================

clear 
close all
clc

ivpfvp = 1;
%% Picard-Chebyshev Convergence Analysis (Example 8: TPBVP)

Nvec = [linspace(10,100,10) linspace(200,1000,9)];
color = {'b*','r*','g*','k*','m*','c*','y*',...
    'bo','ro','go','ko','mo','co','yo',...
    'b^','r^','g^','k^','m^','c^','y^'};

for i = 1:length(Nvec)
    
    N = Nvec(i);
    M = N;
    
    [T2_i,P2_i,T1_i,P1_i,Ta_i,A_i] = clenshaw_curtis_ivpII(N,M);
    [T2_f,P2_f,T1_f,P1_f,Ta_f,A_f] = clenshaw_curtis_fvpII(N,M);
    
    if ivpfvp == 1
        E = eig(T2_f*P2_f*P1_i*A_i);
    elseif ivpfvp == 0;
        E = eig(T2_i*P2_i*P1_f*A_f);
    end
    
    figure(1)
    semilogx(N,E,color{i})
    hold on
    grid on
    
    maxE(i) = max(abs(E));
end

%% Compute tspan for various c and N
cvec = [1e-6 1e-5 1e-4 1e-3 1e-2 0.1 1];
ccolor = {'b*-', 'r*-','g*-','k*-','m*-','c*-','y*-'};
for i = 1:length(cvec)
    c          = cvec(i);
    tspan(:,1) = 2./(cvec(i).*sqrt(abs(maxE)))';
    
    figure(2)
    loglog(Nvec,tspan./60,ccolor{i},'Linewidth',2)
    hold on
end

%% Plots
% N vs max eigenvalue of [T2*PB*P2*P1*A]
figure(1)
grid on
set(gca, 'FontName', 'Helvetica','FontSize',16)
ylabel('|\lambda_{max}(TPBP_2P_1A)|')
xlabel('N')
% Time span for various c and N
figure(2)
grid on
set(gca, 'FontName', 'Helvetica','FontSize',16)
xlabel('N')
ylabel('Time span (mins)')
legend('c = 1e-6','c = 1e-5','c = 1e-4','c = 1e-3','c = 1e-2','c = 0.1','c = 1')

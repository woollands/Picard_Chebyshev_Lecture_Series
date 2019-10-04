% Robyn Woollands 2017
% Texas A&M University - Department of Aerospace Engineering
% File name     : lecture3_example7_ivpII_conv.m
% Description   : Eigenvalue convergence analysis for Lecture 3 of the 
%                 JUNKINS & WOOLLANDS, Picard-Chebyshev Lecture Series.
% Date Written  : March 31, 2017
% Date Modified : March 31, 2017
%================================================================

clear 
close all
clc

%% Picard-Chebyshev Convergence Analysis (Example 7: Second Order IVP)

Nvec = [linspace(10,100,10) linspace(200,1000,9)];
color = {'b*','r*','g*','k*','m*','c*','y*',...
    'bo','ro','go','ko','mo','co','yo',...
    'b^','r^','g^','k^','m^','c^','y^'};

for i = 1:length(Nvec)
    
    N = Nvec(i);
    M = N;

    [T2,P2,T1,P1,Ta,A] = clenshaw_curtis_ivpII(N,M);
    
    E = eig(T2*P2*P1*A);
    
    if N <= 100;
        figure(1)
        hold on
        plot(real(E),imag(E),color{i})
    end
    
   if N > 90; 
        figure(2)
        hold on
        plot(real(E),imag(E),color{i})
    end

    maxE(i) = max(abs(E));
end

%% Compute tspan for various c and N
cvec = [1e-6 1e-5 1e-4 1e-3 1e-2 0.1 1];
ccolor = {'bs-', 'rs-','gs-','ks-','ms-','cs-','ys-'};
for i = 1:length(cvec)
    c          = cvec(i);
    tspan(:,1) = 2./(cvec(i).*sqrt(abs(maxE)))';
    
    figure(4)
    loglog(Nvec,tspan./60,ccolor{i},'Linewidth',2)
    hold on
end

%% Plots
% Root Locus
figure(1)
grid on
axis equal
set(gca, 'FontName', 'Helvetica','FontSize',16)
xlabel('Re(\lambda)')
ylabel('Im(\lambda)')
legend('N=10','N=20','N=30','N=40','N=50',...
    'N=60','N=70','N=80','N=90','N=100')
% Root Locus
figure(2)
grid on
axis equal
set(gca, 'FontName', 'Helvetica','FontSize',16)
xlabel('Re(\lambda)')
ylabel('Im(\lambda)')
legend('N=100','N=200','N=300','N=400','N=500','N=600',...
    'N=700','N=800','N=900','N=1000')
% N vs max eigenvalue of [T2*P2*P1*A]
figure(3)
loglog(Nvec,maxE,'ks-','Linewidth',2)
grid on
set(gca, 'FontName', 'Helvetica','FontSize',16)
xlabel('N')
ylabel('|\lambda_{max}(TP_2P_1A)|')
% Time span for various c and N
figure(4)
grid on
set(gca, 'FontName', 'Helvetica','FontSize',16)
xlabel('N')
ylabel('Time span (mins)')
legend('c = 1e-6','c = 1e-5','c = 1e-4','c = 1e-3','c = 1e-2','c = 0.1','c = 1')

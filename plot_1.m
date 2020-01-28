clear all; clc
r = 10;
s = 10;

alpha = [1:10];
M = 100*alpha;
W = 1000*alpha;

M1=[30,30,30,30,30,30,30,30,30,30];
R1=15;
M2 = M-M1; %for sparsity
R2 = floor(2.9*s*log10(1024/s)); % for sparsity

samp = M.*W;
%scatter(alpha,samp,'filled')
%grid on;
%Fit = polyfit(alpha,samp,2); % x = x data, y = y data, 1 = order of the polynomial i.e a straight line 
%scatter(alpha,samp,80,'filled')
%hold on
%plot(polyval(Fit,alpha),'LineWidth',2.5)
plot(alpha,samp,'-s','MarkerSize',6,'LineWidth',2.5)
hold on

samp1 = R1*M1 + R2.*M2
plot(alpha,samp1,'-o','MarkerSize',6,'LineWidth',2.5)
%Fit = polyfit(alpha,samp1,1)
%scatter(alpha,samp1,80,'filled')
%hold on
%plot(polyval(Fit,alpha),'LineWidth',2.5)
hold on

samp2 = r*((M+W).*log10(M.*W))
plot(alpha,samp2,'-x','MarkerSize',8,'LineWidth',2.5)
%Fit = polyfit(alpha,samp2,1)
%scatter(alpha,samp2,80,'filled')
%hold on
%plot(polyval(Fit,alpha),'LineWidth',2.5)
ylim([-1*10^6 10*10^6])
grid on
%scatter(alpha,samp2,'filled')
ax = gca;
ax.XAxis.FontSize = 11;
ax.YAxis.FontSize = 11;

title({'$$  M=30 \alpha, W=1024 \alpha, S=10, R=10, M_1=20$$','$$\Omega=15, \Delta=60$$'},'fontsize',15,'interpreter','latex')
xlabel( '${\alpha}$','fontsize',24,'interpreter','latex')
ylabel('$T$','fontsize',24,'interpreter','latex')%('Sampling Rate (\Omega M_1 + \Delta M_2)','fontsize',18)

h = legend('Nyquisit Rate','Sparse and Low Rank','Low Rank','Location', 'Best');
h.FontSize = 16;



%a = get(gca,'XTickLabel');
%set(gca,'XTickLabel',a,'FontName','Times','fontsize',10)
%set(gca,'XTickLabel',a, 'FontSize',16)
%set(gca,'FontSize',20)
clc;
close all;
clear

Ks = 33.4*10^9;




rous1 = 2650;
rous2 = 2650;
rous = 2650;
% Kd = 3.18*10^9;


Kd1 =3.18*10^9;
Kd2 =3.18*10^9;

mu1=1.4*10^9;
mu2=1.4*10^9;

fai1 = 0.3;
fai2 = 0.3;


k1=0.5*10^(-13);
k2=0.5*10^(-13);

OO = 1;
F =50;% (fai)^(-1.75);
nj = 8;
% % Kf1= 1*10^6;
% % yita1 = 1*10^(-5);
% % rouf1 = 100;

Kf2= 9.6*10^6;
yita2 = 1.5*10^(-5);
rouf2 = 70;

Kf1= 2.2*10^9;
yita1 = 0.6*10^(-3);
rouf1 = 1000;


% % Kf2= 2.2*10^9;
% % yita2 = 0.6*10^(-3);
% % rouf2 = 1000;

Nx = 400;
Nz = 400;

E = zeros(Nz,Nx);
lamda=zeros(Nz,Nx);
mu = zeros(Nz,Nx);
alfa  = zeros(Nz, Nx);
M_1 = zeros(Nz, Nx);
k0_n = zeros(Nz,Nx);
rouf=zeros(Nz,Nx);



for ix =1:Nx
    for iz = 1:Nz
        
        lamda(iz, ix) = Kd2-2/3*mu2;
        mu(iz,ix)=mu2;
        alfa(iz,ix)=1-Kd2/Ks;
          M_1(iz, ix) = (fai2/Kf2+(alfa(iz,ix)-fai2)/Ks);
        k0_n(iz,ix)= k2/yita2; 
        E(iz, ix) = Kd2+4/3*mu2;
         rouf(iz, ix) =rouf2;
    end
end

for ix =101:300
    for iz =101:300
        lamda(iz, ix) = Kd1-2/3*mu1;
        mu(iz,ix)=mu1;
        alfa(iz,ix)=1-Kd1/Ks;
           M_1(iz, ix) = (fai1/Kf1+(alfa(iz,ix)-fai1)/Ks);       
        k0_n(iz,ix)= k1/yita1; 
         E(iz, ix) = Kd1+4/3*mu1;
          rouf(iz, ix) =rouf1;

    end
end
M_11=M_1;


v =10*10/20/20;



Cmap_min=0.35;
Cmap_max=0.91;
figure;
% cmap=[ 0.25 0.25 0.25;0.81 0.81 0.81];
cmap1=[ Cmap_min Cmap_min Cmap_min;Cmap_max Cmap_max Cmap_max];
imagesc(M_11);
colormap(cmap1)
set(gca, 'units','centimeters','position', [2 1.5 6 6]);
set(gcf,'Units','centimeters','Position',[1 1 13 8]);
axis equal
axis([0 400 0 400]);


dz = 0.001/2/2;
dx = dz;
N_total = Nx*Nz;
V=(Nx-1)*dx*(Nz-1)*dz;
%%   [ 1 t_xx, 2 t_zz, 3 t_xz, 4 P, 5 Us_x, 6 Us_z, 7 Uf_x, 8 Uf_z] 

L1=[-1,1];
L2=[-2,-1,1,2];
%% parpool(4, 'IdleTimeout', 1200)
F_min=-1; F_max=7;
F_inter=6; F_step=1/6;
F_total=F_max-F_min;
kk=0;
for nf=1:F_total
     f_left=F_min+nf-1;
     for ii=1:F_inter
         kk=kk+1;
         ff1(kk)=10^(f_left+(ii-1)*F_step);
     end
end
for nf=1:4
    kk=kk+1;
    ff1(kk)=ff1(kk-1)*10^(2*F_step);
end
nw=length(ff1(:));
rou=v*(rous1*(1-fai1)+fai1*rouf1)+(1-v)*(rous1*(1-fai1)+fai1*rouf2);

alfa1 = 1-Kd1/Ks;
alfa2 = 1-Kd2/Ks;

H1 = Kd1+4/3*mu1+alfa1^2*((fai1/Kf1+(alfa1-fai1)/Ks))^(-1);
H2 = Kd2+4/3*mu2+alfa2^2*((fai2/Kf2+(alfa2-fai2)/Ks))^(-1);
H_hill(1:nw) = ((1-v)/H2+v/H1)^(-1)-4/3*mu1;

%******rewrite 2019-07-17*******%


a_w=1-Kd1/Ks;
Kfw=(v/Kf1+(1-v)/Kf2)^(-1);
H_wood(1:nw)=Kd1+a_w^2*(fai1/Kfw+(a_w-fai1)/Ks)^(-1);


 
 
 
load('KK_1941.mat');
load('KK_1956.mat');



KK=KK_1941;

QQ1 = imag(KK+4/3*mu1)./(real(KK+4/3*mu1));
QQ2 = imag(KK_1956+4/3*mu1)./(real(KK_1956+4/3*mu1));
kk=1;
  figure; 
semilogx(ff1(1:kk:end),QQ1(1:kk:end),'*k','LineWidth',2.5);
hold on;
semilogx(ff1(1:kk:end),QQ2(1:kk:end),'k','LineWidth',2.5);
legend('Biot1941','Biot1956','Location','northwest','fontsize',11,'fontweight','b');
grid on
set(gca,'GridLineStyle','-','GridColor','k','FontName','Times New Roman','FontSize',11,'fontweight','b','LineWidth',1.5)
set(gca, 'units','centimeters','position', [2 1.5 10 6]);
set(gcf,'Units','centimeters','Position',[1 1 13 8]);
 legend('boxoff')
    xlabel('Frequency /Hz','fontsize',12,'fontweight','b','color','k');
    ylabel('Q\^(-1)','fontsize',12,'fontweight','b','color','k');
    axis([ 10^(-1) 10^(8) 0 0.06])
%    print Fig01_b.png -dpng -r800



  kk=1;
figure;
semilogx(ff1(1:kk:end), sqrt((H_hill+4/3*mu1)/rou),':k','LineWidth',2.5);
hold on
semilogx(ff1(1:kk:end), real(sqrt((KK-1/3*mu1+4/3*mu1)/rou)),'*k','LineWidth',2);
hold on
semilogx(ff1(1:kk:end), real(sqrt((KK_1956-1/3*mu1+4/3*mu1)/rou)),'k','LineWidth',2);
hold on
semilogx(ff1(1:kk:end), sqrt((H_wood+4/3*mu1)/rou),'--k','LineWidth',2);
legend('GH-limit','Biot1941','Biot1956','GW-limit','Location','west','fontsize',11,'fontweight','b');
 legend('boxoff')
axis([ 10^(-1) 10^(8) 1600 1750])
set(gca,'GridLineStyle','-','GridColor','k','FontName','Times New Roman','FontSize',11,'fontweight','b','LineWidth',1.5)
grid on
set(gca, 'units','centimeters','position', [2 1.5 10 6]);
set(gcf,'Units','centimeters','Position',[1 1 13 8]);
% % % set(gca,'LooseInset',get(gca,'TightInset'));
% % set(gca, 'Position', get(gca, 'OuterPosition') - ...
% %     get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
    xlabel('Frequency /Hz','fontsize',12,'fontweight','b','color','k');
    ylabel('Velocity /m.s^(^-^1^)','fontsize',12,'fontweight','b','color','k');
% print Fig01_a.png -dpng -r800














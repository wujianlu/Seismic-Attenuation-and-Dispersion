clc;
close all;
clear

Ks = 33.4*10^9;
mu = 1.4*10^9;
rous1 = 2650;
rous2 = 2650;
rous = 2650;
Kd = 3.18*10^9;
fai = 0.3;
k =0.5*10^(-13);
k0=k;

OO = 1;
F =50;% (fai)^(-1.75);
nj = 8;
Kf2= 2.2*10^9;
yita2= 0.6*10^(-3);
rouf2 = 1000;

% % Kf1= 2.2*10^9;
% % yita1 = 0.6*10^(-3);
% % rouf1 = 1000;

Kf1= 9.6*10^6;
yita1 = 1.5*10^(-5);
rouf1 = 70;


dz = 0.0002;
dlayer= 0.1;
n_element = 1;

N_layer = dlayer/dz;
N_total = n_element *dlayer*2/dz;
N_stress = N_total +2;
N_vel = N_total+1;

%%M_stress = zeres(N_stress,1);  %% the peremater of  model  for stress and pressure
%%M_vel = zeros(N_vel,1);        %% the peremater of  model  for velocity


alfa1 = zeros(N_stress,1);
M1 = zeros(N_stress,1);
E1 = zeros(N_stress,1);
E2 = zeros(N_stress,1);
n_k1 = zeros(N_stress,1);
A011 = zeros(N_stress,1);
A021 = zeros(N_stress,1);
Rouf1 = zeros(N_stress,1);

alfa = zeros(N_total,1);
M = zeros(N_total,1);
E = zeros(N_total,1);
E0 = zeros(N_total,1);
n_k = zeros(N_total,1);
A01 = zeros(N_total,1);
A02 = zeros(N_total,1);
Rouf = zeros(N_total,1);
for ki=1:1:n_element
   
    for j = 1:1:N_layer
      l = (ki-1)*N_layer*2+1;
      alfa1(l+j,1) = 1-Kd/Ks;
      M1 (l+j,1) = (fai/Kf1+(alfa1(l+j,1)-fai)/Ks)^(-1);
      n_k1(l+j,1) = yita1 / k0;
      A011(l+j,1) = -rouf1*(1+OO)*F;
      A021(l+j,1) = yita1/k0;
      Rouf1(l+j,1) = rouf1;
      E1(l+j,1) = Kd+4/3*mu;
      E2(l+j,1) = Kd+4/3*mu+alfa1(l+j,1)^2*M1(l+j,1);
      
      s = (ki-1)*N_layer*2+N_layer+1;
     
      alfa1(s+j,1) = 1-Kd/Ks;
      M1 (s+j,1) = (fai/Kf2+(alfa1(s+j,1)-fai)/Ks)^(-1);
      n_k1(s+j,1) = yita2 / k0;
      A011(s+j,1) = -rouf2*(1+OO)*F;
      A021(s+j,1) = yita2/k0;
      Rouf1(s+j,1) = rouf2;     
      E1(s+j,1) = Kd+4/3*mu;
      E2(s+j,1) = Kd+4/3*mu+alfa1(s+j,1)^2*M1(s+j,1);
    end  
    
    
end


for  j = 2:1:(2*n_element*N_layer+1-N_layer/2)
    
     
      E1(j,1) =E1(j+N_layer/2,1);
      E2(j,1) =E2(j+N_layer/2,1);
      alfa1(j,1) =alfa1(j+N_layer/2,1);
      M1(j,1) =M1(j+N_layer/2,1);
      n_k1(j,1) =n_k1(j+N_layer/2,1) ;
      A011(j,1) = A011(j+N_layer/2,1);
      A021(j,1) = A021(j+N_layer/2,1);
      Rouf1(j,1) =Rouf1(j+N_layer/2,1);       
     
end
for  j = (2*n_element*N_layer+2-N_layer/2):1:(2*n_element*N_layer+1)
    
       ki =j-(2*n_element*N_layer-N_layer/2);
    
       E1(j,1) =E1(ki,1);
       E2(j,1) =E2(ki,1);
       alfa1(j,1) =alfa1(ki,1);
       M1(j,1) =M1(ki,1);
       n_k1(j,1) =n_k1(ki,1);
      n_k1(j,1) =n_k1(ki,1) ;
      A011(j,1) = A011(ki,1);
      A021(j,1) = A021(ki,1);
      Rouf1(j,1) =Rouf1(ki,1);        
end





  E1(1,1) =E1(2,1);
   E2(1,1) =E2(2,1);
  alfa1(1,1) =alfa1(2,1);
  M1(1,1) =M1(2,1);
  n_k1(1,1) =n_k1(2,1); 
  A011(1,1) = A011(2,1);
  A021(1,1) = A021(2,1);
  Rouf1(1,1) =Rouf1(2,1);  
  
  
  
  E1(N_stress,1) =E1(N_stress-1,1);
  E2(N_stress,1) =E2(N_stress-1,1); 
  alfa1(N_stress,1) =alfa1(N_stress-1,1);
  M1(N_stress,1) =M1(N_stress-1,1);
  n_k1(N_stress,1) =n_k1(N_stress-1,1);
  A011(N_stress,1) = A011(N_stress-1,1);
  A021(N_stress,1) = A021(N_stress-1,1);
  Rouf1(N_stress,1) =Rouf1(N_stress-1,1);  
    
  
  E0(:,1) =E1(2:(N_total+1),1);
  E(:,1) =E2(2:(N_total+1),1);
  alfa(:,1) =alfa1(2:(N_total+1),1);
  M(:,1) =M1(2:(N_total+1),1);
  n_k(:,1) =n_k1(2:(N_total+1),1);
  A01(:,1) = A011(2:(N_total+1),1);
  A02(:,1) = A021(2:(N_total+1),1);
  Rouf(:,1) =Rouf1(2:(N_total+1),1);  
  
  
  
  
  [EE0,rou,x1,ff1,Vp1,Q1,dt1,dt2]=Biot1941_function(Ks,mu,rous1,rous2,Kd,fai,k,F,Kf2,yita2,rouf2,Kf1,yita1,rouf1,dz,dlayer,n_element,N_layer,N_total,N_stress,N_vel,E0,alfa,M,n_k);
 [EE1,rou2,x2,ff2,Vp2,Q2,dt3, dt4]=Biot1956_function(Ks,mu,rous1,rous2,Kd,fai,k0,F,OO,nj,Kf2,yita2,rouf2,Kf1,yita1,rouf1,dz,dlayer,n_element,N_layer,N_total,N_stress,N_vel,E,alfa,M,n_k,Rouf);

a_w=1-Kd/Ks;
Kfw=(1/2/Kf1+1/2/Kf2)^(-1);
K_wood=Kd+4/3*mu+a_w^2*(fai/Kfw+(a_w-fai)/Ks)^(-1);

Kh1=Kd+4/3*mu+a_w^2*(fai/Kf1+(a_w-fai)/Ks)^(-1);
Kh2=Kd+4/3*mu+a_w^2*(fai/Kf2+(a_w-fai)/Ks)^(-1);
K_hill=(1/2/Kh1+1/2/Kh2)^(-1);

N_f=length(ff1);
Vp_hill(1:N_f)=sqrt((K_hill)/rou);
Vp_wood(1:N_f)=sqrt(K_wood/rou);
  kk=1;
figure;
semilogx(ff1(1:kk:end), Vp1(1:kk:end),'k','LineWidth',2);
hold on
semilogx(ff1(1:kk:end), sqrt(abs(real(EE0(1:kk:end))/rou)),'o','Color',[0.8500 0.3250 0.0980],'LineWidth',1.5)
hold on
semilogx(ff1(1:kk:end), sqrt(real(EE1(1:kk:end))/rou),'Color',[0 0.4470 0.7410],'LineWidth',2)
hold on
semilogx(ff1(1:kk:end), Vp_wood(1:kk:end),'--k','LineWidth',2);
hold on
semilogx(ff1(1:kk:end), Vp_hill(1:kk:end),':k','LineWidth',2.5);
legend('White1975','Biot1941','Biot1956','GW-limit','GH-limit','Location','west','fontsize',11,'fontweight','b');
 legend('boxoff')
axis([ 10^(-2) 10^(7) 1550 1870])
set(gca,'GridLineStyle','-','GridColor','k','FontName','Times New Roman','FontSize',11,'fontweight','b','LineWidth',1.5)

grid on
set(gca, 'units','centimeters','position', [2 1.5 15 6]);
set(gcf,'Units','centimeters','Position',[1 1 18 8]);
% % % set(gca,'LooseInset',get(gca,'TightInset'));
% % set(gca, 'Position', get(gca, 'OuterPosition') - ...
% %     get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
    xlabel('Frequency /Hz','fontsize',12,'fontweight','b','color','k');
    ylabel('Velocity /m.s^(^-^1^)','fontsize',12,'fontweight','b','color','k');
% print Fig01_a.png -dpng -r800

% grid on
% %     
% % set (gcf, 'visible','off')
% % 
% % print(gcf,'-r500', 'Vp0.eps', '-depsc');


kk=1;
  figure; 
semilogx(ff1(1:kk:end),Q1(1:kk:end),'k','LineWidth',2.5);
 hold on
  semilogx(ff1(1:kk:end),imag(EE0(1:kk:end))./real(EE0(1:kk:end)),'o','Color',[0.8500 0.3250 0.0980],'LineWidth',2)
hold on
semilogx(ff1(1:kk:end),imag(EE1(1:kk:end))./real(EE1(1:kk:end)),'-','Color',[0 0.4470 0.7410],'LineWidth',2)
legend('White1975','Biot1941','Biot1956','Location','northwest','fontsize',11,'fontweight','b');
grid on
set(gca,'GridLineStyle','-','GridColor','k','FontName','Times New Roman','FontSize',11,'fontweight','b','LineWidth',1.5)
set(gca, 'units','centimeters','position', [2 1.5 15 6]);
set(gcf,'Units','centimeters','Position',[1 1 18 8]);
 legend('boxoff')
    xlabel('Frequency /Hz','fontsize',12,'fontweight','b','color','k');
    ylabel('Q\^(-1)','fontsize',12,'fontweight','b','color','k');
    axis([ 10^(-2) 10^(7) 0 0.14])
%    print Fig01_b.png -dpng -r800

   
  Z= (1:1:N_total)*dz;
  X = ff1;   

  %% displacement
 figure;
 colormap('hot')
 
 h=surf(X,Z,abs(x1(1+N_total:2*N_total,:)));
  set(gca,'XScale','log')
  set(h,'edgecolor','none');
  shading interp;
set(gca,'GridLineStyle','-','GridColor','k','FontName','Times New Roman','FontSize',11,'fontweight','b','LineWidth',1.5)
 set(gcf,'Units','centimeters','Position',[1 1 15 8]);%����ͼƬ��СΪ13cm��10cm
     xlabel('Frequency /Hz','fontsize',12,'fontweight','b','color','k');
    ylabel('Thickness /m','fontsize',12,'fontweight','b','color','k');
    colorbar('FontName','Times New Roman','FontSize',11,'FontWeight','b');
      set(gca,'ydir','reverse')
          axis([ 10^(-2) 10^(7) 0 0.2])
% %     CC=colorbar;%%
% %      CC.Label.FontSize=11;
% %      CC.Label.FontWeight='b';
    view(0,90)
% print Fig02_a.png -dpng -r800


%% Fluid pressure

   
figure;
 colormap('hot')
 
 h=surf(X,Z,abs(x1(1+3*N_total:4*N_total,:)));
  set(gca,'XScale','log')
  set(h,'edgecolor','none');
  shading interp;
set(gca,'GridLineStyle','-','GridColor','k','FontName','Times New Roman','FontSize',11,'fontweight','b','LineWidth',1.5)
 set(gcf,'Units','centimeters','Position',[1 1 15 8]);%����ͼƬ��СΪ13cm��10cm
     xlabel('Frequency /Hz','fontsize',12,'fontweight','b','color','k');
    ylabel('Thickness /m','fontsize',12,'fontweight','b','color','k');
    colorbar('FontName','Times New Roman','FontSize',11,'FontWeight','b');
        axis([ 10^(-2) 10^(7) 0 0.2])
         set(gca,'ydir','reverse')
    view(0,90)
% print Fig02_b.png -dpng -r800
   
  
  
  
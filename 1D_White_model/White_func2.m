
%% ���ǵ����ͺ�ˮ��

% % % clc;
% % % clear;
% % % close;
% % % S1=0.3;S2=0.3;   % ��д��Ϊ��϶��%
% % % Ks1=33.4*10^9;Ks2=33.4*10^9;    %GPa%
% % % Kf1=2.2*10^9;Kf2=6*10^5;  %GPa%
% % % d1=0.1;d2=0.1;    %m%
% % % n1=0.003; %% ˮ��ճ����
% % % n2=1.5*10^(-5);  %GPa.s%
% % % %k1=1000/9.81*10^(-14);k2=100/9.81*10^(-14); 
% % % k1=1*10^(-13);k2=1*10^(-13);
% % % 
% % % raus1=2650;rauf1=1000;
% % % rauf2=80;
% % % %%  Gurevich parameater
% % % Km=6*10^9;
% % % Kb=10*10^9;
% % % b_R=0.001;
% % % faic=0.0002;
% % % mu=6*10^9;

function[Vp,Q,ff,dt1, dt2 ]=White_func2(S1,S2,Ks1,Ks2,Kf1,Kf2,d1,d2,n1,n2,k1,k2,raus1,rauf1,rauf2,Kb,mu1,n)
% n=50;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%��1��Ĳ���%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%% The following is  the besic parameters  of the first phase %%%%%%%%%
%% 



nw= 180;
for nn = 1:1:nw
    
    aa=(-2+80/280*(nn-1)*0.2);
    %ff(f)=log(10^aa)/log(10);
    ff(nn)=10^(aa);

     w=ff(nn)*pi*2;   
  w1=nn;   

K2=Kb;
u2=mu1;
K1=K2;
u1=u2;


%%


a1=1-K1/Ks1;a2=1-K2/Ks2;




ka1=(S1/Kf1+(1-S1)/Ks1-K1/Ks1^2)^(-1);ka2=(S2/Kf2+(1-S2)/Ks2-K2/Ks2^2)^(-1);
M1=K1+4/3*u1+a1^2*ka1;M2=K2+4/3*u2+a2^2*ka2;
h1=d1/(d1+d2);h2=d2/(d1+d2);
C0=(h1/M1+h2/M2)^(-1);  
Md1=K1+4/3*u1;Md2=K2+4/3*u2;
R1=a1*ka1/M1;R2=a2*ka2/M2;
ke1=ka1*Md1/M1;ke2=ka2*Md2/M2;
Z1=(n1*ke1./(i.*w*k1)).^(1/2).*coth((i*w*n1/(k1*ke1)).^(1/2)*d1);Z2=(n2*ke2./(i*w*k2)).^(1/2).*coth((i*w*n2/(k2*ke2)).^(1/2)*d2);
%%  ��ȷ���ʽ
C=C0./(1+C0*(R1-R2)^2./((Z1+Z2)*i.*w.*(d1+d2)));
C1(w1,1)=C;
Ci=imag(C);
Cr=real(C);
Q(w1,1)=Ci./Cr;
Vp(w1,1)=real(sqrt(C/(raus1*(1-S1)+rauf1*S2*h1+rauf2*S2*h2)));
f(w1,1)=real(C)-K1;
rau=raus1*(1-S1)+rauf1*S2*h1+rauf2*S2*h2;
%%
Sg=d1/(d1+d2);

K_w=(Sg/Kf1+(1-Sg)/Kf2)^(-1);
K_BGW(w1,1)=K1+a1^2*((a1-S1)/Ks1+S1/K_w)^(-1)+4/3*u1;
K_BG1=K1+a1^2*((a1-S1)/Ks1+S1/Kf1)^(-1)+4/3*u1;
K_BG2=K1+a1^2*((a1-S1)/Ks1+S1/Kf2)^(-1)+4/3*u1;
K_BGH(w1,1)=(Sg/K_BG1+(1-Sg)/K_BG2)^(-1);

end


c1 = 9/8;
c2 = 1/24;
dx =0.01;
rou = (raus1*(1-S1)+rauf1*S2*h1+rauf2*S2*h2);
Vp1 = sqrt((K_BG1/rou));

Vp2 = sqrt((K_BG1/rou));
dt1=dx/(c1+c2)/Vp1;
dt2=dx/(c1+c2)/Vp2;




% %  semilogx(ff,Vp);
% % hold on;
% % semilogx(K_BGW,'r')
% % hold on;
% % semilogx(K_BGH,'--r')




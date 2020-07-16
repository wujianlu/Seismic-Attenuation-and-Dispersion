%% Using rubino(2009)'s boundary condition ,Solving the White model in frequncy domine
function [EE0,rou,x1,ff1,Vp1,Q1,dt1, dt2]=Biot1941_function(Ks,mu,rous1,rous2,Kd,fai,k,F,Kf2,yita2,rouf2,Kf1,yita1,rouf1,dz,dlayer,n_element,N_layer,N_total,N_stress,N_vel,E,alfa,M,n_k)





% % % % Ks = 40*10^9;
% % % % mu = 3*10^9;
% % % % rous1 = 2700;
% % % % rous2 = 2700;
% % % % rous = 2700;
% % % % Kd = 4*10^9;
% % % % fai = 0.2;
% % % % k =1*10^(-13);
% % % % 

% % % % F =100;% (fai)^(-1.75);
% % % % 
% % % % Kf2= 2.3*10^9;
% % % % yita2 = 3*10^(-3);
% % % % rouf2 = 1000;
% % % % 
% % % % Kf1= 0.022*10^9;
% % % % yita1 = 1*10^(-5);
% % % % rouf1 = 140;
% % % % 
% % % % dz = 0.001;
% % % % dlayer= 0.1;
% % % % n_element = 1;

% % % % N_layer = dlayer/dz;
% % % % N_total = n_element *dlayer*2/dz;
% % % % N_stress = N_total +2;
% % % % N_vel = N_total+1;
% % % % 
% % % % %%M_stress = zeres(N_stress,1);  %% the peremater of  model  for stress and pressure
% % % % %%M_vel = zeros(N_vel,1);        %% the peremater of  model  for velocity
% % % % 
% % % % 
% % % % 
% % % % alfa1 = zeros(N_stress,1);
% % % % M1 = zeros(N_stress,1);
% % % % E1 = zeros(N_stress,1);
% % % % n_k1 = zeros(N_stress,1);
% % % % alfa = zeros(N_total,1);
% % % % M = zeros(N_total,1);
% % % % E = zeros(N_total,1);
% % % % n_k = zeros(N_total,1);
% % % % 
% % % % for ki=1:1:n_element
% % % %    
% % % %     for j = 1:1:N_layer
% % % %       l = (ki-1)*N_layer*2+1;
% % % %     
% % % %       E1(l+j,1) = Kd+4/3*mu;
% % % %       alfa1(l+j,1) = 1-Kd/Ks;
% % % %       M1 (l+j,1) = (fai/Kf1+(alfa1(l+j,1)-fai)/Ks)^(-1);
% % % %       n_k1(l+j,1) = yita1 / k;
% % % %       
% % % %      
% % % %       
% % % %       s = (ki-1)*N_layer*2+N_layer+1;
% % % %       E1(s+j,1) = Kd+4/3*mu;
% % % %       alfa1(s+j,1) = 1-Kd/Ks;
% % % %       M1 (s+j,1) = (fai/Kf2+(alfa1(s+j,1)-fai)/Ks)^(-1);
% % % %       n_k1(s+j,1) = yita2 / k;
% % % %       
% % % %       
% % % %     end  
% % % %     
% % % %     
% % % % end
% % % % 
% % % % 
% % % % for  j = 2:1:(2*n_element*N_layer+1-N_layer/2)
% % % %     
% % % %      
% % % %       E1(j,1) =E1(j+N_layer/2,1);
% % % %       alfa1(j,1) =alfa1(j+N_layer/2,1);
% % % %       M1(j,1) =M1(j+N_layer/2,1);
% % % %       n_k1(j,1) =n_k1(j+N_layer/2,1) ;
% % % %       
% % % %      
% % % % end
% % % % for  j = (2*n_element*N_layer+2-N_layer/2):1:(2*n_element*N_layer+1)
% % % %     
% % % %        ki =j-(2*n_element*N_layer-N_layer/2);
% % % %     
% % % %        E1(j,1) =E1(ki,1);
% % % %        alfa1(j,1) =alfa1(ki,1);
% % % %        M1(j,1) =M1(ki,1);
% % % %        n_k1(j,1) =n_k1(ki,1);
% % % %        
% % % % end
% % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % %   E1(1,1) =E1(2,1);
% % % %   alfa1(1,1) =alfa1(2,1);
% % % %   M1(1,1) =M1(2,1);
% % % %   n_k1(1,1) =n_k1(2,1);
% % % %   
% % % %   E1(N_stress,1) =E1(N_stress-1,1);
% % % %   alfa1(N_stress,1) =alfa1(N_stress-1,1);
% % % %   M1(N_stress,1) =M1(N_stress-1,1);
% % % %   n_k1(N_stress,1) =n_k1(N_stress-1,1);
% % % %   
% % % %   
% % % %   E(:,1) =E1(2:(N_total+1),1);
% % % %   alfa(:,1) =alfa1(2:(N_total+1),1);
% % % %   M(:,1) =M1(2:(N_total+1),1);
% % % %   n_k(:,1) =n_k1(2:(N_total+1),1);
  
  
  
 tao=1;


c1= 9/8;
c2 = 1/24;
N_N = (N_total+N_vel*2)*(N_total+N_vel*2);
nw= 180;
for nn = 1:1:nw
    
    aa=(-2+80/280*(nn-1)*0.2);
    %ff(f)=log(10^aa)/log(10);
    ff1(nn)=10^(aa);
    f0=ff1(nn);
    iw =i*2*pi*f0;
    index = 0;
        
    %% equation 111111111  div t = 0
     for iz = 1:N_total
         for m =1:2
                if((iz+2*m-3)>0&&(iz+2*m-3)<(N_total+1))
                        index = index+1;
                        sparse_i(1,index)=iz;
                        sparse_j(1,index)=iz+2*m-3;   
                       switch (  m  )
                           case {1}
                                   sparse_v(1,index)=-1/dz/2;
                           case {2}
                                   sparse_v(1,index)=1/dz/2;
                       end
                end
         end
     
     end
    %%  equation 2222222   t+ap-E*Dz*Us=0
    for iz =1:1:N_total
         index = index+1;
         sparse_i(1,index)=iz+N_total;
         sparse_j(1,index)=iz;  
         sparse_v(1,index) =1;
         
         index = index+1;
         sparse_i(1,index)=iz+N_total;
         sparse_j(1,index)=iz+N_total;  
         sparse_v(1,index) =alfa(iz,1);
      
     
         
         for m = 1:1:2
          
             if((iz+2*m-3)>0&&(iz+2*m-3)<(N_total+1))   
                  index = index+1;
                  sparse_i(1,index)=iz+N_total;
                  sparse_j(1,index)=iz+2*m-3+2*N_total;  
                  switch (  m  )
                      case {1}
                           sparse_v(1,index)=E(iz,1)/dz/2;
                       case {2}
                           sparse_v(1,index)=-E(iz,1)/dz/2;
                  end            
             end          
         end
         
         
         
     end
    %% equation 3 3    Dz*Uf+a*Dz*Us+P/M=0
    for iz =1:1:N_total
     for m = 1:1:2          
             if((iz+2*m-3)>0&&(iz+2*m-3)<(N_total+1))   
                  index = index+1;
                  sparse_i(1,index)=iz+2*N_total;
                  sparse_j(1,index)=iz+2*m-3+3*N_total;  
                  switch (  m  )
                      case {1}
                           sparse_v(1,index)=-1/dz/2;
                       case {2}
                           sparse_v(1,index)=1/dz/2;
                  end            
             end          
          end
          for m = 1:1:2          
             if((iz+2*m-3)>0&&(iz+2*m-3)<(N_total+1))   
                  index = index+1;
                  sparse_i(1,index)=iz+2*N_total;
                  sparse_j(1,index)=iz+2*m-3+2*N_total;  
                  switch (  m  )
                      case {1}
                           sparse_v(1,index)=-alfa(iz,1)/dz/2;
                       case {2}
                           sparse_v(1,index)=alfa(iz,1)/dz/2;
                  end            
          end      
     end
          
         index = index+1;
         sparse_i(1,index)=iz+2*N_total;
         sparse_j(1,index)=iz+N_total; 
         sparse_v(1,index) =1/M(iz,1);   
         
    
    end
  
    %% equation 444444  iw*Uf+k/n*Dz*p=0  
    
     for iz =1:1:N_total
         
         index = index+1;
         sparse_i(1,index)=iz+3*N_total;
         sparse_j(1,index)=iz+3*N_total;  
         sparse_v(1,index) =iw;
         
          for m = 1:1:2          
             if((iz+2*m-3)>0&&(iz+2*m-3)<(N_total+1))   
                  index = index+1;
                  sparse_i(1,index)=iz+3*N_total;
                  sparse_j(1,index)=iz+2*m-3+1*N_total;  
                  switch (  m  )
                      case {1}
                           sparse_v(1,index)=-1/n_k(iz,1)/dz/2;
                       case {2}
                           sparse_v(1,index)=1/n_k(iz,1)/dz/2;
                  end            
             end          
          end          
    
     end
  
 %%  calculate AA    
     
     AA = sparse( sparse_i,sparse_j, sparse_v);
    
     b = zeros(4*N_total,1);
     
% % % %   %% add boundary   stress =t0   (t0 =1); 
% % % %      AA(:,(1+N_total)) = 0;
% % % %     AA((1+N_total),:) = 0;
% % % %     AA(:,(2*N_total)) = 0;
% % % %     AA((2*N_total),:) = 0;
% % % %     
% % % %     AA((1+N_total),(1+N_total)) =1;
% % % %     AA(2*N_total,2*N_total) = 1;
    
    %% Add    boundary bottom   Uf  =0;

    AA(:,(1+3*N_total)) = 0;
    AA((1+3*N_total),:) = 0;
    AA(1+3*N_total,1+3*N_total) = 1;
    AA(:,(4*N_total)) = 0;
    AA((4*N_total),:) = 0;
    AA(4*N_total,4*N_total) = 1;
% % % %     %% Add    boundary bottom   Us  =0;
% % % % 
% % % %  
% % % %     AA(:,(3*N_total)) = 0;
% % % %     AA((3*N_total),:) = 0;
% % % %     AA(3*N_total,3*N_total) = 1;    
    
    
    
    b(1,1) = 1;
   b(N_total,1) = -1;
    
    
    n =0;
    n1= 0;
       
    
    
    
    
    %%
    dd= AA\b;
   x1(:,nn)=dd;
    
    
    
    
    
    
    
    
    
    
end

ee = 0;
for nn1 =1:nw
    st =0;
    
    ss0 = 0;
    for nn = (1+ee):(1-ee)
          st =st+x1(nn,nn1);
          ss0 =ss0+x1(nn+1+2*N_total,nn1)/2/dz; 
  
    end
    
    for nn = (2+ee):(N_total-1)
          st =st+x1(nn,nn1);
          ss0 =ss0+(x1(nn+1+2*N_total,nn1)-x1(nn-1+2*N_total,nn1))/2/dz; 
  
    end   
    for nn = (N_total):(N_total)
          st =st+x1(nn,nn1);
          ss0 =ss0-x1(nn-1+2*N_total,nn1)/2/dz; 
  
    end   
    stress(nn1,1) = st;
    strain(nn1,1) = ss0;
    
    
end
rou=0.5*(rous1*(1-fai)+fai*rouf1)+0.5*(rous1*(1-fai)+fai*rouf2);

EE0= stress./strain;
[Vp1,Q1,ff1,dt1, dt2]=White_func2(fai,fai,Ks,Ks,Kf1,Kf2,dlayer/2,dlayer/2,yita1,yita2,k,k,rous1,rouf1,rouf2,Kd,mu,nw);

detu = (-x1(1+2*N_total,:)+x1(3*N_total,:))/(n_element *dlayer*2);
% EE0 = 4*dz./detu;

% % % % figure;
% % % % semilogx(sqrt(real(EE0)/rou),'or')
% % % % hold on
% % % % semilogx(ff1, Vp1,'k','LineWidth',2);
% % % % legend('��ֵ��','White1975');
% % % % set(gcf,'Units','centimeters','Position',[10 10 10 6]);%����ͼƬ��СΪ13cm��10cm
% % % % grid on
% % % % 
% % % % 
% % % % figure;
% % % % semilogx(imag(EE0)./real(EE0),'or')
% % % % hold on;
% % % % semilogx(ff1,Q1,'k','LineWidth',2);
% % % % legend('��ֵ��','White1975');
% % % % set(gcf,'Units','centimeters','Position',[10 10 10 6]);%����ͼƬ��СΪ13cm��10cm
% % % % grid on




























clc;
close all;
clear

Ks = 33.4*10^9;

rous1 = 2650;
rous2 = 2650;
rous = 2650;


Kd1 =3.18*10^9;
Kd2 =3.18*10^9;

mu1=1.4*10^9;
mu2=1.4*10^9;

fai1 = 0.3;
fai2 = 0.3;


k1=0.5*10^(-13);
k2=0.5*10^(-13);

OO = 1;
F =50;
nj = 8;

Kf2= 9.6*10^6;
yita2 = 1.5*10^(-5);
rouf2 = 70;

Kf1= 2.2*10^9;
yita1 = 0.6*10^(-3);
rouf1 = 1000;

Nx = 120;
Nz = 120;

E = zeros(Nz,Nx);
lamda=zeros(Nz,Nx);
mu = zeros(Nz,Nx);
alfa  = zeros(Nz, Nx);
M_1 = zeros(Nz, Nx);
k0_n = zeros(Nz,Nx);
rouf=zeros(Nz,Nx);



for ix =1:Nx
    for iz = 1:Nz        
        mu(iz,ix)=mu2;
        alfa(iz,ix)=1-Kd2/Ks;
          M_1(iz, ix) = (fai2/Kf2+(alfa(iz,ix)-fai2)/Ks);
        k0_n(iz,ix)= k2/yita2; 
        E(iz, ix) = Kd2+4/3*mu2+alfa(iz,ix)*alfa(iz,ix)/M_1(iz,ix);
        rouf(iz, ix) =rouf2;
        lamda(iz, ix) = Kd2-2/3*mu2+alfa(iz,ix)*alfa(iz,ix)/M_1(iz,ix);
    end
end

for ix =31:90
    for iz = 31:90
        
        mu(iz,ix)=mu1;
        alfa(iz,ix)=1-Kd1/Ks;
           M_1(iz, ix) = (fai1/Kf1+(alfa(iz,ix)-fai1)/Ks);       
        k0_n(iz,ix)= k1/yita1; 
         E(iz, ix) = Kd1+4/3*mu1+alfa(iz,ix)*alfa(iz,ix)/M_1(iz,ix);
         lamda(iz, ix) = Kd1-2/3*mu1+alfa(iz,ix)*alfa(iz,ix)/M_1(iz,ix);
        rouf(iz, ix) =rouf1;
    end
end
dz = 0.005/6;
dx = dz;
N_total = Nx*Nz;
V=(Nx-1)*dx*(Nz-1)*dz;

v =10*10/20/20;

index=0;
    %% equation 111111111  Dx*t_xx+Dz*t_xz+w^2*rous*Us_x+w^2*rou_f*Uf_x  = 0    
    for ix=1:Nx
        for iz = 1:Nz
            for m=1:2
                 if((iz+2*m-3)>0&&(iz+2*m-3)<(Nz+1))
                        index = index+1;
                        sparse_i(1,index)=(iz-1)*Nx+ix;
                        sparse_j(1,index)=(iz+2*m-3-1)*Nx+ix+2*N_total;   
                       switch (  m  )
                           case {1}
                                   sparse_v(1,index)=-1/dz/2;
                           case {2}
                                   sparse_v(1,index)=1/dz/2;
                       end                      
                       
                  end            
                  if((ix+2*m-3)>0&&(ix+2*m-3)<(Nx+1))
                        index = index+1;
                        sparse_i(1,index)=(iz-1)*Nx+ix;
                        sparse_j(1,index)=(iz-1)*Nx+ix+2*m-3;   
                       switch (  m  )
                           case {1}
                                   sparse_v(1,index)=-1/dx/2;
                           case {2}
                                   sparse_v(1,index)=1/dx/2;
                       end                      
                  end                    
            end
        end
    end

      %% equation 222222222222   Dx*t_xz+Dz*t_zz +w^2*rous*Us_z+w^2*rou_f*Uf_z  = 0 
     for ix=1:Nx
        for iz = 1:Nz
            for m=1:2
                 if((iz+2*m-3)>0&&(iz+2*m-3)<(Nz+1))
                        index = index+1;
                        sparse_i(1,index)=(iz-1)*Nx+ix+N_total;
                        sparse_j(1,index)=(iz+2*m-3-1)*Nx+ix+N_total;   
                       switch (  m  )
                           case {1}
                                   sparse_v(1,index)=-1/dz/2;
                           case {2}
                                   sparse_v(1,index)=1/dz/2;
                       end                      
                  end            
                  if((ix+2*m-3)>0&&(ix+2*m-3)<(Nx+1))
                        index = index+1;
                        sparse_i(1,index)=(iz-1)*Nx+ix+N_total;
                        sparse_j(1,index)=(iz-1)*Nx+ix+2*m-3+2*N_total;   
                       switch (  m  )
                           case {1}
                                   sparse_v(1,index)=-1/dx/2;
                           case {2}
                                   sparse_v(1,index)=1/dx/2;
                       end                      
                  end                         
            end
        end
     end

    
    
%% equation 333333333333   t_xx - E* Dx*Us_x - lamda*Dz*Us_z  -a*M*Dx*Uf_x-a*M*Dz*Uf_z=0      
   for ix=1:Nx
        for iz = 1:Nz
            index = index+1;
            sparse_i(1,index)=(iz-1)*Nx+ix+2*N_total;
            sparse_j(1,index)=(iz-1)*Nx+ix;  
            sparse_v(1,index) =1;
% %             index = index+1;
% %             sparse_i(1,index)=(iz-1)*Nx+ix+2*N_total;
% %             sparse_j(1,index)=(iz-1)*Nx+ix+3*N_total;  
% %             sparse_v(1,index) =alfa(iz,ix);           
            
            for m=1:2  
                 if((iz+2*m-3)>0&&(iz+2*m-3)<(Nz+1))
                        index = index+1;
                        sparse_i(1,index)=(iz-1)*Nx+ix+2*N_total;
                        sparse_j(1,index)=(iz+2*m-3-1)*Nx+ix+5*N_total;   
                       switch (  m  )
                           case {1}
                                   sparse_v(1,index)=lamda(iz,ix)/dz/2;
                           case {2}
                                   sparse_v(1,index)=-lamda(iz,ix)/dz/2;
                       end                      
                  end  
                  if((ix+2*m-3)>0&&(ix+2*m-3)<(Nx+1))
                        index = index+1;
                        sparse_i(1,index)=(iz-1)*Nx+ix+2*N_total;
                        sparse_j(1,index)=(iz-1)*Nx+ix+2*m-3+4*N_total;   
                       switch (  m  )
                           case {1}
                                   sparse_v(1,index)=E(iz,ix)/dx/2;
                           case {2}
                                   sparse_v(1,index)=-E(iz,ix)/dx/2;
                       end                      
                  end                  

                 if((iz+2*m-3)>0&&(iz+2*m-3)<(Nz+1))
                        index = index+1;
                        sparse_i(1,index)=(iz-1)*Nx+ix+2*N_total;
                        sparse_j(1,index)=(iz+2*m-3-1)*Nx+ix+7*N_total;   
                       switch (  m  )
                           case {1}
                                   sparse_v(1,index)=alfa(iz,ix)/dz/2/M_1(iz,ix);
                           case {2}
                                   sparse_v(1,index)=-alfa(iz,ix)/dz/2/M_1(iz,ix);
                       end                      
                  end  
                  if((ix+2*m-3)>0&&(ix+2*m-3)<(Nx+1))
                        index = index+1;
                        sparse_i(1,index)=(iz-1)*Nx+ix+2*N_total;
                        sparse_j(1,index)=(iz-1)*Nx+ix+2*m-3+6*N_total;   
                       switch (  m  )
                           case {1}
                                   sparse_v(1,index)=alfa(iz,ix)/dz/2/M_1(iz,ix);
                           case {2}
                                   sparse_v(1,index)=-alfa(iz,ix)/dz/2/M_1(iz,ix);
                       end                      
                  end                  

            end
            
        end
    end          
   %% equation 444444444444   t_zz - E* Dz*Us_z  - lamda*Dx*Us_x  -a*M*Dx*Uf_x-a*M*Dz*Uf_z=0      
   for ix=1:Nx
        for iz = 1:Nz
            index = index+1;
            sparse_i(1,index)=(iz-1)*Nx+ix+3*N_total;
            sparse_j(1,index)=(iz-1)*Nx+ix+N_total;  
            sparse_v(1,index) =1;
% %             index = index+1;
% %             sparse_i(1,index)=(iz-1)*Nx+ix+3*N_total;
% %             sparse_j(1,index)=(iz-1)*Nx+ix+3*N_total;  
% %             sparse_v(1,index) =alfa(iz,ix);           
            
            for m=1:2          
                 if((iz+2*m-3)>0&&(iz+2*m-3)<(Nz+1))
                        index = index+1;
                        sparse_i(1,index)=(iz-1)*Nx+ix+3*N_total;
                        sparse_j(1,index)=(iz+2*m-3-1)*Nx+ix+5*N_total;   
                       switch (  m  )
                           case {1}
                                   sparse_v(1,index)=E(iz,ix)/dz/2;
                           case {2}
                                   sparse_v(1,index)=-E(iz,ix)/dz/2;
                       end                      
                 end    
                  if((ix+2*m-3)>0&&(ix+2*m-3)<(Nx+1))
                        index = index+1;
                        sparse_i(1,index)=(iz-1)*Nx+ix+3*N_total;
                        sparse_j(1,index)=(iz-1)*Nx+ix+2*m-3+4*N_total;   
                       switch (  m  )
                           case {1}
                                   sparse_v(1,index)=lamda(iz,ix)/dx/2;
                           case {2}
                                   sparse_v(1,index)=-lamda(iz,ix)/dx/2;
                       end                      
                  end
                  if((iz+2*m-3)>0&&(iz+2*m-3)<(Nz+1))
                        index = index+1;
                        sparse_i(1,index)=(iz-1)*Nx+ix+3*N_total;
                        sparse_j(1,index)=(iz+2*m-3-1)*Nx+ix+7*N_total;   
                       switch (  m  )
                           case {1}
                                   sparse_v(1,index)=alfa(iz,ix)/dz/2/M_1(iz,ix);
                           case {2}
                                   sparse_v(1,index)=-alfa(iz,ix)/dz/2/M_1(iz,ix);
                       end                      
                  end  
                  if((ix+2*m-3)>0&&(ix+2*m-3)<(Nx+1))
                        index = index+1;
                        sparse_i(1,index)=(iz-1)*Nx+ix+3*N_total;
                        sparse_j(1,index)=(iz-1)*Nx+ix+2*m-3+6*N_total;   
                       switch (  m  )
                           case {1}
                                   sparse_v(1,index)=alfa(iz,ix)/dx/2/M_1(iz,ix);
                           case {2}
                                   sparse_v(1,index)=-alfa(iz,ix)/dx/2/M_1(iz,ix);
                       end                      
                  end                     
              
            end
            
        end
    end           
   %% equation 555555555555555  t_xz-mu*(Dx*Us_z + Dz*Us_x)=0      
    for ix=1:Nx
        for iz = 1:Nz
            index = index+1;
            sparse_i(1,index)=(iz-1)*Nx+ix+4*N_total;
            sparse_j(1,index)=(iz-1)*Nx+ix+2*N_total;  
            sparse_v(1,index) =1;         
            
            for m=1:2          
                 if((iz+2*m-3)>0&&(iz+2*m-3)<(Nz+1))
                        index = index+1;
                        sparse_i(1,index)=(iz-1)*Nx+ix+4*N_total;
                        sparse_j(1,index)=(iz+2*m-3-1)*Nx+ix+4*N_total;   
                       switch (  m  )
                           case {1}
                                   sparse_v(1,index)=mu(iz,ix)/dz/2;
                           case {2}
                                   sparse_v(1,index)=-mu(iz,ix)/dz/2;
                       end                      
                  end            
                  if((ix+2*m-3)>0&&(ix+2*m-3)<(Nx+1))
                        index = index+1;
                        sparse_i(1,index)=(iz-1)*Nx+ix+4*N_total;
                        sparse_j(1,index)=(iz-1)*Nx+ix+2*m-3+5*N_total;   
                       switch (  m  )
                           case {1}
                                   sparse_v(1,index)=mu(iz,ix)/dx/2;
                           case {2}
                                   sparse_v(1,index)=-mu(iz,ix)/dx/2;
                       end                      
                  end                                       
            end            
        end
    end        
    %% equation 666666666666666  Dx*Uf_x + Dz*Uf_z +a*(Dx*Us_x + Dz*Uz_z) +M_1*P=0    
    for ix=1:Nx
        for iz = 1:Nz
            index = index+1;
            sparse_i(1,index)=(iz-1)*Nx+ix+5*N_total;
            sparse_j(1,index)=(iz-1)*Nx+ix+3*N_total;  
            sparse_v(1,index) =M_1(iz,ix);         
            
            for m=1:2          
                 if((iz+2*m-3)>0&&(iz+2*m-3)<(Nz+1))
                        index = index+1;
                        sparse_i(1,index)=(iz-1)*Nx+ix+5*N_total;
                        sparse_j(1,index)=(iz+2*m-3-1)*Nx+ix+7*N_total;   
                       switch (  m  )
                           case {1}
                                   sparse_v(1,index)=-1/dz/2;
                           case {2}
                                   sparse_v(1,index)=1/dz/2;
                       end                      
                  end            
                  if((ix+2*m-3)>0&&(ix+2*m-3)<(Nx+1))
                        index = index+1;
                        sparse_i(1,index)=(iz-1)*Nx+ix+5*N_total;
                        sparse_j(1,index)=(iz-1)*Nx+ix+2*m-3+6*N_total;   
                       switch (  m  )
                           case {1}
                                    sparse_v(1,index)=-1/dx/2;
                           case {2}
                                     sparse_v(1,index)=1/dx/2;
                       end                      
                  end                    
                    
            end
             for m=1:2          
                 if((iz+2*m-3)>0&&(iz+2*m-3)<(Nz+1))
                        index = index+1;
                        sparse_i(1,index)=(iz-1)*Nx+ix+5*N_total;
                        sparse_j(1,index)=(iz+2*m-3-1)*Nx+ix+5*N_total;   
                       switch (  m  )
                           case {1}
                                   sparse_v(1,index)=-alfa(iz,ix)/dz/2;
                           case {2}
                                   sparse_v(1,index)=alfa(iz,ix)/dz/2;
                       end                      
                  end            
                  if((ix+2*m-3)>0&&(ix+2*m-3)<(Nx+1))
                        index = index+1;
                        sparse_i(1,index)=(iz-1)*Nx+ix+5*N_total;
                        sparse_j(1,index)=(iz-1)*Nx+ix+2*m-3+4*N_total;   
                       switch (  m  )
                           case {1}
                                    sparse_v(1,index)=-alfa(iz,ix)/dx/2;
                           case {2}
                                     sparse_v(1,index)=alfa(iz,ix)/dx/2;
                       end                      
                  end                                        
            end           
        end
    end     
     %% equation 777777777777 i*w*Uf_x+k_n*Dx*P +w^2*rou_f*Us_x=0    

     for ix=1:Nx
        for iz = 1:Nz   
            for m=1:2                
                  if((ix+2*m-3)>0&&(ix+2*m-3)<(Nx+1))
                        index = index+1;
                        sparse_i(1,index)=(iz-1)*Nx+ix+6*N_total;
                        sparse_j(1,index)=(iz-1)*Nx+ix+2*m-3+3*N_total;   
                       switch (  m  )
                           case {1}
                                    sparse_v(1,index)=-k0_n(iz,ix)/dx/2;
                           case {2}
                                    sparse_v(1,index)=k0_n(iz,ix)/dx/2;
                       end                      
                  end                                      
            end          
        end
    end       
       %% equation 888888888888 i*w*Uf_z+k_n*Dz*P +w^2*rou_f*Us_z=0   
     for ix=1:Nx
        for iz = 1:Nz    
            for m=1:2                 
                 if((iz+2*m-3)>0&&(iz+2*m-3)<(Nz+1))
                        index = index+1;
                        sparse_i(1,index)=(iz-1)*Nx+ix+7*N_total;
                        sparse_j(1,index)=(iz+2*m-3-1)*Nx+ix+3*N_total;   
                       switch (  m  )
                           case {1}
                                    sparse_v(1,index)=-k0_n(iz,ix)/dz/2;
                           case {2}
                                    sparse_v(1,index)=k0_n(iz,ix)/dz/2;
                       end                      
                  end                                       
            end         
        end
     end   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
     for ix=1:Nx
        for iz = 1:Nz           
            index = index+1;
            sparse_i(1,index)=(iz-1)*Nx+ix;
            sparse_j(1,index)=(iz-1)*Nx+ix+4*N_total; 
            sparse_v(1,index) = rous;  
            index = index+1;
            sparse_i(1,index)=(iz-1)*Nx+ix;
            sparse_j(1,index)=(iz-1)*Nx+ix+6*N_total; 
            sparse_v(1,index) = rouf(iz,ix);  
            
        end
     end     
      for ix=1:Nx
        for iz = 1:Nz    
             index = index+1;
            sparse_i(1,index)=(iz-1)*Nx+ix+N_total;
            sparse_j(1,index)=(iz-1)*Nx+ix+5*N_total; 
            sparse_v(1,index) = rous;  
            index = index+1;
            sparse_i(1,index)=(iz-1)*Nx+ix+N_total;
            sparse_j(1,index)=(iz-1)*Nx+ix+7*N_total; 
            sparse_v(1,index) = rouf(iz,ix);          
        end
     end  
     for ix=1:Nx
        for iz = 1:Nz
             index = index+1;
            sparse_i(1,index)=(iz-1)*Nx+ix+6*N_total;
            sparse_j(1,index)=(iz-1)*Nx+ix+4*N_total;  
            sparse_v(1,index) =rouf(iz,ix)*k0_n(iz,ix);    
             index = index+1;
            sparse_i(1,index)=(iz-1)*Nx+ix+7*N_total;
            sparse_j(1,index)=(iz-1)*Nx+ix+5*N_total;  
            sparse_v(1,index) =rouf(iz,ix)*k0_n(iz,ix); 
        end
     end    
     for ix=1:Nx
        for iz = 1:Nz
            index = index+1;
            sparse_i(1,index)=(iz-1)*Nx+ix+6*N_total;
            sparse_j(1,index)=(iz-1)*Nx+ix+6*N_total;  
            sparse_v(1,index) =1;   
            index = index+1;
            sparse_i(1,index)=(iz-1)*Nx+ix+7*N_total;
            sparse_j(1,index)=(iz-1)*Nx+ix+7*N_total;  
            sparse_v(1,index) =1;         
        end
     end
N_length=index;

%%   [ 1 t_xx, 2 t_zz, 3 t_xz, 4 P, 5 Us_x, 6 Us_z, 7 Uf_x, 8 Uf_z] 
% parpool(4, 'IdleTimeout', 1200)
F_min=-1; F_max=7;
F_inter=6; F_step=1/6;
F_total=F_max-F_min;
kk=0;
for nf=1:F_total
     f_left=F_min+nf-1;
     for ii=1:F_inter
         kk=kk+1;
         ff2(kk)=10^(f_left+(ii-1)*F_step);
     end
end
aa=ff2(end);
ff1(1)=aa*10^(2*F_step);
kk=1;
for nf=1:4
    kk=kk+1;
    ff1(kk)=ff1(kk-1)*10^(2*F_step);
end
 nw=length(ff1(:));    
for nn=1:nw
    

    f0=ff1(nn)
    iw =i*2*pi*f0;
    w2=2*pi*f0*2*pi*f0;
%     index = 0;   
    sparse_vv=zeros(1,N_length);
    sparse_vv=sparse_v;
    sparse_vv(1,((N_length-6*Nx*Nz+1):(N_length-2*Nx*Nz)))=sparse_v(1,((N_length-6*Nx*Nz+1):(N_length-2*Nx*Nz)))*w2;
    sparse_vv(1,((N_length-2*Nx*Nz+1):(N_length)))=sparse_v(1,((N_length-2*Nx*Nz+1):(N_length)))*iw;   
     AA = sparse( sparse_i,sparse_j, sparse_vv);
    
     b = zeros(8*N_total,1);    

     
  %% Us_z=0  Uf_z=0  P=0  on the right and left side
   ix_L = 1;
   ix_R = Nx;

   for iz =1:1:Nz
       

            AA(ix_L+(iz-1)*Nx+6*N_total,:) = 0;
            AA(:,ix_L+(iz-1)*Nx+6*N_total) = 0;                   
            AA(ix_R+(iz-1)*Nx+6*N_total,:) = 0;
            AA(:,ix_R+(iz-1)*Nx+6*N_total) = 0; 
            AA(ix_L+(iz-1)*Nx+6*N_total,ix_L+(iz-1)*Nx+6*N_total) = 1;   
            AA(ix_R+(iz-1)*Nx+6*N_total,ix_R+(iz-1)*Nx+6*N_total) = 1;   

            
   end

   %%   Us_x=0;  Uf_x = 0  P=0 on the top and bottom side 
   iz_T = 1;
   iz_B = Nz;   
   for ix =1:1:Nx
       
      
      AA(:, (iz_T-1)*Nx+ix+7*N_total) = 0; 
      AA((iz_T-1)*Nx+ix+7*N_total,:) = 0; 
      AA(:, (iz_B-1)*Nx+ix+7*N_total) = 0; 
      AA((iz_B-1)*Nx+ix+7*N_total,:) = 0; 
      AA((iz_T-1)*Nx+ix+7*N_total,(iz_T-1)*Nx+ix+7*N_total) = 1; 
      AA((iz_B-1)*Nx+ix+7*N_total,(iz_B-1)*Nx+ix+7*N_total) = 1;    
      
   end
     

   
   
   
   for ix=1:Nx
       b(N_total+ix,1) = 1; 
       b(N_total+(iz_B-1)*Nx+ix,1) = -1; 
   end
   
 
       
   for iz=1:Nz
       b((iz-1)*Nx+ix_L,1) =1; 
       b((iz-1)*Nx+ix_R,1) = -1; 
   end   


   
% %    for iz=1:Nz
% %        b(N_total+(iz-1)*Nx+ix_L,1) =1; 
% %        b(N_total+(iz-1)*Nx+ix_R,1) = -1; 
% %    end   
   
   
    
    
    dd= AA\b;
   x1(:,1)=dd;
   temp=0;
   for ix=1:Nx
      temp = temp +  x1(ix+5*N_total,1);
   end  
  Us_z_top_average (nn,1) = temp/Nx;
   temp=0; 
   for ix=1:Nx
      temp = temp +  x1(ix+(Nz-1)*Nx+5*N_total,1);
   end  
  Us_z_bottom_average (nn,1) = temp/Nx;
    
    temp=0;
   for iz=1:Nz
      temp = temp +  x1(1+(iz-1)*Nx+4*N_total,1);
   end  
  Us_x_left_average (nn,1) = temp/Nx;
    temp=0;
   for iz=1:Nz
      temp = temp +  x1(Nx+(iz-1)*Nx+4*N_total,1);
   end  
  Us_x_right_average (nn,1) = temp/Nx;
  
  det_V = (V-((Nx-1)*dx-(Us_x_right_average(nn,1)-Us_x_left_average(nn,1))*((Nz-1)*dz-(Us_z_bottom_average(nn,1)-Us_z_top_average(nn,1)))));
  
  
  K(nn,1) = 10^3/det_V*V;
  Q(nn,1) = imag(K(nn,1))/real(K(nn,1));
   
   

  for ix = 1: Nx
      for iz = 1: Nz
         Txx(iz,ix)  = x1((iz-1)*Nx+ix,1);
         Tzz(iz,ix)  = x1((iz-1)*Nx+ix+N_total,1);
         Txz(iz,ix)  = x1((iz-1)*Nx+ix+2*N_total,1);
         P(iz,ix)    = x1((iz-1)*Nx+ix+3*N_total,1); 
         Us_x(iz,ix) = x1((iz-1)*Nx+ix+4*N_total,1);
         Us_z(iz,ix) = x1((iz-1)*Nx+ix+5*N_total,1);
         Uf_x(iz,ix) = x1((iz-1)*Nx+ix+6*N_total,1);
         Uf_z(iz,ix) = x1((iz-1)*Nx+ix+7*N_total,1);
      end
  end
  strain_xx = zeros(Nz,Nx);
  strain_zz = zeros(Nz,Nx);
for ix=1:1
    for iz = 1:Nz
       strain_xx(iz,ix)= (Us_x(iz,ix+1))/2/dx;        
    end
end
for ix=2:(Nx-1)
    for iz = 1:Nz
       strain_xx(iz,ix)= (Us_x(iz,ix+1)-Us_x(iz,ix-1))/2/dx;        
    end
end  
 for ix=Nx:Nx
    for iz = 1:Nz
       strain_xx(iz,ix)= (-Us_x(iz,ix-1))/2/dx;        
    end
 end 

for ix=1:Nx
    for iz = 1:1
       strain_zz(iz,ix)= (Us_z(iz+1,ix))/2/dz;        
    end
end
 for ix=1:Nx
    for iz = 2:(Nz-1)
       strain_zz(iz,ix)= (Us_z(iz+1,ix)-Us_z(iz-1,ix))/2/dz;        
    end
end
  for ix=1:Nx
    for iz = Nz:Nz
       strain_zz(iz,ix)= (-Us_z(iz-1,ix))/2/dz;        
    end
  end
 
temp1 =0;
temp2 = 0;
temp3 =0;
temp4 = 0;

% % % % % % % % % % % % % % % % % % % % %   for ix=1:Nx
% % % % % % % % % % % % % % % % % % % % %     for iz = 1:Nz
% % % % % % % % % % % % % % % % % % % % %         temp1 = temp1+Tzz(iz,ix)+Txx(iz,ix);
% % % % % % % % % % % % % % % % % % % % %         temp2 = temp2+strain_zz(iz,ix)+strain_xx(iz,ix);   
% % % % % % % % % % % % % % % % % % % % %         temp3 = temp3+Tzz(iz,ix)-Txx(iz,ix);
% % % % % % % % % % % % % % % % % % % % %         temp4 = temp4+strain_zz(iz,ix)-strain_xx(iz,ix);           
% % % % % % % % % % % % % % % % % % % % %     end
% % % % % % % % % % % % % % % % % % % % %   end
  for ix=1:Nx
    for iz = 1:Nz
% %         temp1 = temp1+Tzz(iz,ix);
% %         temp2 = temp2+strain_zz(iz,ix);   
        temp1 = temp1+Tzz(iz,ix)+Txx(iz,ix);
        temp2 = temp2+strain_zz(iz,ix)+strain_xx(iz,ix);              
        temp3 = temp3+Tzz(iz,ix)-Txx(iz,ix);
        temp4 = temp4+strain_zz(iz,ix)-strain_xx(iz,ix);           
    end
  end  
  
%   temp1 = temp1+1000*2*2*Nx;
  KK_1956(nn,1) = (temp1)/temp2/2;%%-1/3*1/2*temp3/temp4; %% Masson 2007

  KK_shear(nn,1) = 1/2*temp3/temp4;  %% Quaital 2011
 

end
save('KK_1956.mat','KK_1956');

KK=KK_1956;
figure;
semilogx(ff1,real(K(:,1)))
figure;
semilogx(ff1,Q)

figure;
imagesc(abs(Txx));
title('Txx');
figure;
imagesc(abs(Tzz));
title('Tzz');
figure;
imagesc(abs(Txz));
title('Txz');
figure;
imagesc(abs(P));
title('P');
figure;
imagesc(abs(Us_x));
title('Us_x');

figure;
imagesc(abs(Us_z));
title('Us_z');


figure;
imagesc(abs(Uf_x));
title('Uf_x');


figure;
imagesc(abs(Uf_z));
title('Uf_z');



QQ1 = imag(KK+4/3*mu1)./(real(KK+4/3*mu1));
QQs = imag(KK_shear)./(real(KK_shear));

figure;

semilogx(ff1, QQ1,'LineWidth',2)
legend('Numberical Results');
set(gcf,'Units','centimeters','Position',[1 1 16 12]);%����ͼƬ��СΪ13cm��10cm
set(gca,'GridLineStyle','-','GridColor','k','LineWidth',1.5)
    xlabel('Frequency /Hz','fontsize',12,'fontweight','b','color','k');
    ylabel('Q\^(-1)','fontsize',12,'fontweight','b','color','k');
    axis([ 10^(-2) 10^(8) 0 0.18])
grid on

figure;

semilogx(ff1, QQs,'LineWidth',2)
legend('Numberical Results');
set(gcf,'Units','centimeters','Position',[10 10 16 8]);%����ͼƬ��СΪ13cm��10cm
set(gca,'GridLineStyle','-','GridColor','k','LineWidth',1.5)
    xlabel('Frequency /Hz','fontsize',12,'fontweight','b','color','k');
    ylabel('Q\^(-1)','fontsize',12,'fontweight','b','color','k');
%     axis([ 10^(-2) 10^(5) 0 0.012])
grid on



% rou_eff = ((1-fai1)*rous+fai1*(rouf1*v+rouf2*(1-v)));
rou_eff=v*(rous1*(1-fai1)+fai1*rouf1)+(1-v)*(rous1*(1-fai1)+fai1*rouf2);


figure
semilogx(ff1,sqrt((real(KK-1/3*mu1+4/3*mu1))/rou_eff),'LineWidth',2);
legend('Numberical Results');
set(gcf,'Units','centimeters','Position',[1 1 16 8]);%����ͼƬ��СΪ13cm��10cm
    xlabel('Frequency /Hz','fontsize',12,'fontweight','b','color','k');
    ylabel('Velocity /m.s\^(-1)','fontsize',12,'fontweight','b','color','k');
    set(gca,'GridLineStyle','-','GridColor','k','LineWidth',1.5)
%     axis([10^(-2) 10^(5) 3230 3300])
% set(gca,'XTick',0:20:100);
% set(gca,'YTick',2500:200:3600);
hold on
xL=xlim;yL=ylim;
plot(xL,[yL(2),yL(2)],'k',[xL(2),xL(2)],[yL(1),yL(2)],'k')
box off
axis([xL yL]) 
grid on

figure
semilogx(ff1,sqrt((real(KK_shear))/rou_eff),'LineWidth',2);
legend('Numberical Results');
set(gcf,'Units','centimeters','Position',[10 10 16 8]);%����ͼƬ��СΪ13cm��10cm
    xlabel('Frequency /Hz','fontsize',12,'fontweight','b','color','k');
    ylabel('Velocity /m.s\^(-1)','fontsize',12,'fontweight','b','color','k');
    set(gca,'GridLineStyle','-','GridColor','k','LineWidth',1.5)
%     axis([10^(-2) 10^(5) 3230 3300])
% set(gca,'XTick',0:20:100);
% set(gca,'YTick',2500:200:3600);
hold on
xL=xlim;yL=ylim;
plot(xL,[yL(2),yL(2)],'k',[xL(2),xL(2)],[yL(1),yL(2)],'k')
box off
axis([xL yL]) 
grid on;




alfa1 = 1-Kd1/Ks;
alfa2 = 1-Kd2/Ks;

H1 = Kd1+4/3*mu1+alfa1^2*((fai1/Kf1+(alfa1-fai1)/Ks))^(-1);
H2 = Kd2+4/3*mu2+alfa2^2*((fai2/Kf2+(alfa2-fai2)/Ks))^(-1);



temp=0;
% H_wood(1:nw) = (1-v)*H2+v*H1-temp;
H_hill(1:nw) = ((1-v)/H2+v/H1)^(-1)-4/3*mu1;

%******rewrite 2019-07-17*******%

Kfw=(v/Kf1+(1-v)/Kf2)^(-1);

a_1=1-Kd1/Ks;
a_2=1-Kd2/Ks;
a_w=(v/a_1+(1-v)/a_2)^(-1);
Kdw=(v/Kd1+(1-v)/Kd2)^(-1);
M_1=(fai1/Kf1+(a_1-fai1)/Ks)^(-1);
M_2=(fai2/Kf2+(a_2-fai2)/Ks)^(-1);
Mw=(v/M_1+(1-v)/M_2)^(-1);
muw=(v/mu1+(1-v)/mu2)^(-1);
H_wood(1:nw)=Kdw+4/3*mu1+a_w^2*Mw-temp;

a_w=1-Kd1/Ks;
Kfw=(v/Kf1+(1-v)/Kf2)^(-1);
H_wood(1:nw)=Kd1+a_w^2*(fai1/Kfw+(a_w-fai1)/Ks)^(-1);




% % % % a_w=1-Kd/Ks;
% % % % Kfw=(1/3/Kf1+1/3/Kf3+1/3/Kf3)^(-1);
% % % % K_wood=Kd+4/3*mu+a_w^2*(fai/Kfw+(a_w-fai)/Ks)^(-1);
% % % % 
% % % % Kh1=Kd+4/3*mu+a_w^2*(fai/Kf1+(a_w-fai)/Ks)^(-1);
% % % % Kh2=Kd+4/3*mu+a_w^2*(fai/Kf2+(a_w-fai)/Ks)^(-1);
% % % % Kh3=Kd+4/3*mu+a_w^2*(fai/Kf3+(a_w-fai)/Ks)^(-1);
% % % % K_hill=(1/3/Kh1+1/3/Kh2+1/3/Kh3)^(-1);



figure;
semilogx(ff1,real(KK-1/3*mu1+4/3*mu1));
hold on;
semilogx(ff1,(H_wood+4/3*mu1),'r');
hold on;
semilogx(ff1,(H_hill+4/3*mu1),'k');
grid on;

Vp=3000;
Vs=1400;


temp=4/3*mu1;
figure;
semilogx(ff1,sqrt(real(KK-1/3*mu1+4/3*mu1)/rou_eff));
hold on;
semilogx(ff1,sqrt((H_wood+4/3*mu1)/rou_eff),'r');
hold on;
semilogx(ff1,sqrt((H_hill+4/3*mu1)/rou_eff),'k');
grid on;

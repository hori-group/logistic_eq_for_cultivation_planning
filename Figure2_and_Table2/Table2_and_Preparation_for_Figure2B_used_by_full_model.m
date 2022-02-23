%%
%% Identification program for MR model for Figure2 (B) and Table2 
%%

clear;
close all;


%% excel file
   filename = '../_OD_excel_data/OD_data.xlsx';
 
   sheet='data'; 
   xlRange='e4:CE20';
   
    d_in=readmatrix(filename,'Sheet',sheet,'Range',xlRange);
     
    time_data=d_in(1,:);
    od_data=d_in(2:17,:) %OD data 


%% initial values
 ODN=30.10*1e6;
 J= 0.1529;
 K=J/ODN %extinction coefficient
 
 dep = 3.7; %depth D (unit is cm)
 vol = 200; %Volume V (unit is mL) 

 OD0=0.025 %initial OD
 N0 = OD0*ODN*vol; %initial OD fixed

 
%% incident light flux into culture
ein=[96.8, 184.4, 386.7, 1034]; %incident light flux density ein (unit is \mu E s^{-1} m^{-2})
light_area=0.002826 %area (unit is m^2)
Ein=ein*light_area %incident light flux into culture (unit is \mu E s^{-1}) 

time_conversion=3600 %unit is hour s^{-1}
    
Ein=Ein*time_conversion % unit is \mu E hour^{-1}

C0=[1/8,1/4,1/2,1]


%% MR model simulation paramter (see Table2)
 
 mu=0.194;
 lambda=1900*(1e-9);
 Gamma=9.9932e+10;
 Gamma_OD= Gamma/(ODN*vol)

 
%% Identification for Table2 and LOOCV for full model
    
    LOOCV_sign=zeros(4,4);
    %       C(0)=1/8 C(0)=1/4 C(0)=1/2 C(0)=1
    %L(-3)
    %L(-2)
    %L(-1)
    %L(0)
    
    %---------Not LOOCV for Table2---------------------%
    comment='Idenfication for Table2'
   
    %------------LOOCV at Ein=0.274------------------%  
    
    if 0 %LOOCV for Ein=0.274 and C0=1/8
    LOOCV_sign(1,1)=1;
    od_data(1,:)=NaN;
    comment='LOOCV for Ein=0.274 and C0=1/8'
    end
    
     if 0 %LOOCV for Ein=0.274 and C0=1/4
    LOOCV_sign(1,2)=1;
    od_data(2,:)=NaN;
    comment='LOOCV for Ein=0.274 and C0=1/4'
    end
    
     if 0 %LOOCV for Ein=0.274 and C0=1/2
    LOOCV_sign(1,3)=1;
    od_data(3,:)=NaN;
    comment='LOOCV for Ein=0.274 and C0=1/2'
     end
    
      if 0 %LOOCV for Ein=0.274 and C0=1
    LOOCV_sign(1,4)=1;
    od_data(4,:)=NaN;
    comment='LOOCV for Ein=0.274 and C0=1'
      end
      
      
    %------------LOOCV at Ein=0.521------------------%
   
     if 0 %LOOCV for Ein=0.521 and C0=1/8
    LOOCV_sign(2,1)=1;
    od_data(5,:)=NaN;
    comment='LOOCV for Ein=0.521 and C0=1/8'
     end
    
    if 0 %LOOCV for Ein=0.521 and C0=1/4
    LOOCV_sign(2,2)=1;
    od_data(6,:)=NaN;
    comment='LOOCV for Ein=0.521 and C0=1/4'
    end
    
    if 0 %LOOCV for Ein=0.521 and C0=1/2
    LOOCV_sign(2,3)=1;
    od_data(7,:)=NaN;
    comment='LOOCV for Ein=0.521 and C0=1/2'
    end
    
    if 0 %LOOCV for Ein=0.521 and C0=1
    LOOCV_sign(2,4)=1;
    od_data(8,:)=NaN;
     comment='LOOCV for Ein=0.521 and C0=1'
     end
 
    %------------LOOCV at Ein=1.09------------------%
    if 0 %LOOCV for Ein=1.09 and C0=1/8
    LOOCV_sign(3,1)=1
     od_data(9,:)=NaN;
    comment='LOOCV for Ein=1.09 and C0=1/8'
    end
    
    if 0 %LOOCV for Ein=1.09 and C0=1/4
    LOOCV_sign(3,2)=1
     od_data(10,:)=NaN;
    comment='LOOCV for Ein=1.09 and C0=1/4'
    end
    
    if 0 %LOOCV for Ein=1.09 and C0=1/2
    LOOCV_sign(3,3)=1
     od_data(11,:)=NaN;
    comment='LOOCV for Ein=1.09 and C0=1/2'
    end
    
    if 0 %LOOCV for Ein=1.09 and C0=1
    LOOCV_sign(3,4)=1
     od_data(12,:)=NaN;
    comment='LOOCV for Ein=1.09 and C0=1'
    end

%-----------L(0)--------------    
    
    if 0 %LOOCV for Ein=2.92 and C0=1/8
    LOOCV_sign(4,1)=1
     od_data(13,:)=NaN;
    comment='LOOCV for Ein=2.92 and C0=1/8'
   end
    
    if 0 %LOOCV for Ein=2.92 and C0=1/4
    LOOCV_sign(4,2)=1
     od_data(14,:)=NaN;
    comment='LOOCV for Ein=2.92 and C0=1/4'
    end
    
    if 0 %LOOCV for Ein=2.92 and C0=1/2
    LOOCV_sign(4,3)=1
     od_data(15,:)=NaN;
    comment='LOOCV for Ein=2.92 and C0=1/2'
    end
    
    if 0 %LOOCV for Ein=2.92 and C0=1
    LOOCV_sign(4,4)=1
     od_data(16,:)=NaN;
    comment='LOOCV for Ein=2.92 and C0=1'
    end

    
%% Weight coefficient for MSE
length(od_data)

num_of_nan_data=sum(isnan(od_data),2) %summation for nan data
   
num_of_data=length(od_data)-num_of_nan_data %summation for data in each condition

sum_of_num_of_data=sum(num_of_data,'all'); 
   
weight_for_mse=num_of_data/sum_of_num_of_data %Weight coefficient

    
%%
%% Main code
%%  
 simtime = 0:1:1500;
 
alpha=[7:0.1:10]*1e-12;
xi=[0.005:0.001:0.05];
    
   
  MSE=zeros(length(alpha),length(xi));
  se_for_E0_and_C0=zeros(length(Ein),length(C0));
  MSE_for_E0_and_C0=zeros(length(Ein),length(C0));
  OD_differ=zeros(length(time_data),1);
  se=zeros(length(time_data),1);
  OD=zeros(length(simtime),length(Ein));
  
  ind=0;
  
 
  
  for c=1:length(xi)
    
         for a=1:length(alpha)
            for S=1:length(Ein)
                
                for  M=1:length(C0)
                 
                        
                     ind=ind+1;
                     
                      if LOOCV_sign(S,M) %LOCCV_sign
                         
                          se_for_E0_and_C0(S,M)=0;
                         
                      else
                         
                          y0=[N0,C0(M)];    
                          [hour, y] = ode23s(@(t,y) Full_model(t,y,alpha(a),mu,lambda,Ein(S),dep,vol,K,Gamma,xi(c),ODN),...
                                simtime, y0);

                           OD(:,M)=y(:,1)./(ODN*vol);
                           
                                                       
                            for j=1:num_of_data(ind) 
                           
                                 z= find(hour == time_data(j)); 
                                 
                                 OD_differ(j)= OD(z,M)-od_data(ind,j) ;
                                 se(j)=(OD_differ(j))^2;

                                
                            end %j end
                                 
                                se_for_E0_and_C0(S,M)= sum(se,'all'); %SE
                                OD_differ=zeros(length(time_data),1);
                                se=zeros(length(time_data),1);
                                
                                                                                  
                      end %LOOCV sign end
                    
                      MSE_for_E0_and_C0(S,M)= weight_for_mse(ind)*se_for_E0_and_C0(S,M); 
                      
                      
                 end %media end
            end  %E0 end
             

           MSE(a,c)= sum(MSE_for_E0_and_C0,'all');
           ind=0;

         end %alpha end

  end %xi end
    
  
  MIN_MSE=min(MSE,[],'all')
  [alpha_index,xi_index]=find(MSE==MIN_MSE)
  
  Identified_alpha=alpha(alpha_index)
  Identified_xi=xi(xi_index)
   
  
 
  
 %% Identified parameter should be as below: 
% 
%   Table2
%     alpha=8.7*1e-12;
%     xi=0.012;
% 
%   Identified parameter matrix by LOOCV
%        row: light fulx into culture, column: initial media concentration
%                   C0=1/8 C0=1/4 C0=1/2 C0=1
%        Ein=0.274 
%        Ein=0.521
%        Ein=1.09
%        Ein=2.92
%            alpha=[8.7, 8.7, 8.7, 8.7;
%                   8.7, 8.7, 8.5, 8.7;
%                   8.7, 8.8, 8.7, 8.7;
%                   8.7, 8.7, 8.6, 8.7]*1e-12;
%    
%
%              xi=[0.015,0.012,0.013,0.012;
%                  0.014,0.013,0.007,0.012;
%                  0.013,0.011,0.01,0.014;
%                  0.012,0.012,0.011,0.01];
%                

  comment
 
%% function that returns dy/dt

function dy=Full_model(t,y,alpha,mu,lambda,Ein,dep,vol,K,Gamma,xi,ODN)
 
  dy(1,:)= mu.*(y(2)/(xi+y(2)))*(lightpercell(Ein,y(1),dep,vol,K,ODN)/(lambda + lightpercell(Ein,y(1),dep,vol,K,ODN)))*(1-y(1)/Gamma)*y(1);

  dy(2,:)= -alpha*dy(1,:);
   
  
end


%% function lightpercell
function LPC=lightpercell(L,N,dep,vol,KK,ODN)

        cell_conc=@(N) N/vol;
        LPC=L*(1 - 10.^(-KK.*cell_conc(N)*dep))./N;

end



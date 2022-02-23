%%
%% Identification program for MR model for Figure2 (C) and Table2 
%%

clear;
close all;


%% excel file
   filename = '../_OD_excel_data/OD_data.xlsx';
 
   sheet='data'; 
   xlRange='e4:CE20';
   
    d_in=readmatrix(filename,'Sheet',sheet,'Range',xlRange);
     
    time_data=d_in(1,:);
    od_data_for_MRmodel=d_in(5:4:17,:) %OD data with C0=1
 

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


%% Identification for Table2 and LOOCV for MR model
   
    LOOCV_sign=zeros(1,4);
    
    %---------Not LOOCV for Table2---------------------%
    comment='Idenfication for Table2'
    
    %----------------LOOCV-----------------------------%
    if 0 %LOOCV for Ein=0.274
    LOOCV_sign(1,1)=1;
    od_data_for_MRmodel(1,:)=NaN;
    comment='LOOCV for Ein=0.274'
    end
    
    if 0 %LOOCV for Ein=0.521
    LOOCV_sign(1,2)=1;
    od_data_for_MRmodel(2,:)=NaN;
     comment='LOOCV for Ein=0.521'
    end
    
    if 0 %LOOCV for Ein=1.09
    LOOCV_sign(1,3)=1;
     od_data_for_MRmodel(3,:)=NaN;
    comment='LOOCV for Ein=1.09'
    end
    
    if 0 %LOOCV for Ein=2.92
    LOOCV_sign(1,4)=1;
    od_data_for_MRmodel(4,:)=NaN;
    comment='LOOCV for Ein=2.92'
     end
    
    
%% Corrected the data used for MR model


index1=30; %306 hour
index2=61; %690 hour

od_data_for_MRmodel(3,index2+1:length(time_data))=NaN; %OD data after 690 hour was altered NaN
od_data_for_MRmodel(4,index1+1:length(time_data))=NaN; %OD data after 306 hour was altered NaN

nun_index=isnan(od_data_for_MRmodel); 
sum_of_nun=sum(nun_index,2); %Sum of NaN data

index=length(time_data); %sum of number of data
od_data_index=[index-sum_of_nun(1),index-sum_of_nun(2),...
    index-sum_of_nun(3),index-sum_of_nun(4)] %fixed number of index 


%% Weight coefficient for MSE

number_of_data=sum(od_data_index) 
    
weight_for_mse=od_data_index/number_of_data %Weight coefficient

%% Main code
            
   if  LOOCV_sign(1,4)~=1
       mu=[0.18:0.001:0.25];
       lambda=[1700:50:2500]*(1e-9);
       Gamma_OD=[15:0.1:19];
       Gamma=Gamma_OD*(ODN*vol);
   end

           
    if  LOOCV_sign(1,4)==1
        mu=[0.2:0.001:0.35];
        lambda=[2000:50:4000]*(1e-9);
        Gamma_OD=[15:0.1:18];
        Gamma=Gamma_OD*(ODN*vol);
    end
            
simtime = 0:1:1500;

  for mu_ind=1:length(mu)
       
    for lambda_ind=1:length(lambda)
        
        for Gamma_ind=1:length(Gamma)
            
for  S=1:length(Ein)
      
    if weight_for_mse(S)==0
        
        MSE_for_each_E0(S)=0;
        
    else
       [hour, y] = ode45(@(t,y) MRmodel(t,y,mu(mu_ind),lambda(lambda_ind),Ein(S),dep,vol,K,Gamma(Gamma_ind),ODN),simtime, N0);
       %y is biomass N
        OD=y(:,1)/(ODN*vol);

                 for j=1:od_data_index(S) 

                             k= find(hour == time_data(j));

                             OD_differ(j)= OD(k)-od_data_for_MRmodel(S,j) ;
                            mse(j)=(OD_differ(j))^2;

                 end
                 
                 se_for_each_E0(S)= sum(mse,'all'); %SE
                 MSE_for_each_E0(S)= weight_for_mse(S)*se_for_each_E0(S); %MSE

                OD_differ=[];
                 mse=[];     
    end % if end
    
 end %slit end
       
        MSE(lambda_ind,Gamma_ind,mu_ind)= sum(MSE_for_each_E0,'all');
                MSE_for_each_E0=[];

        end %Gamma end 
    end %lambda end 
  end %mu end 
   
   %MSE
  MIN_MSE=min(MSE,[],'all')
  [lambda_index,linear_index]=find(MSE==MIN_MSE)
  sz=[length(Gamma),length(mu)];
  [Gamma_index,mu_index] = ind2sub(sz,linear_index)
  
  
   Identified_mu=mu(mu_index)
   
   Identified_lambda=lambda(lambda_index)
   
   Identified_Gamma=Gamma(Gamma_index)
   Identified_Gamma_OD=Identified_Gamma/(ODN*vol)
   
%% end   

comment
 
%% Identified parameter should be as below: 
% 
%   Table2
%     mu=0.194;
%     lambda=1900*(1e-9);
%    Gamma=9.9932e+10;
%    Gamma_OD=16.6
% 
%  
%   LOOCV for Ein=0.274
%      mu=0.187;
%      lambda=1850*(1e-9);
%      Gamma= 1.0174e+11;
%      Gamma_OD=16.9
%
%   LOOCV for Ein=0.521
%      mu=0.182;
%      lambda=1750*(1e-9);
%      Gamma=  9.8728e+10;
%      Gamma_OD=16.4
%    
%   LOOCV for Ein=1.09
%      mu=0.181;
%      lambda=1750*(1e-9);
%      Gamma=  9.7524e+10;
%      Gamma_OD=16.2
%
%   LOOCV for Ein=2.92
%      mu=0.349;
%      lambda=3550*(1e-9);
%      Gamma=  1.0174e+11;
%      Gamma_OD=16.9
%  

 
%% function that returns dy/dt
 
function dy = MRmodel(t,y,mu,lambda,E0,dep,vol,K,Gamma,ODN)

dy = mu.*(lightpercell(E0,y,dep,vol,K,ODN)/(lambda + lightpercell(E0,y,dep,vol,K,ODN)))*(1-y/Gamma)*y;
        
end


%% function lightpercell
function LPC=lightpercell(L,N,dep,vol,K,ODN)

        cell_conc=@(N) N/vol;
        LPC=L*(1 - 10.^(-K.*cell_conc(N)*dep))./N;

end


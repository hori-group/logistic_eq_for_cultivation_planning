%%
%% Figure2 (C)  
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


   
%%  
%% main code for Figure2C
%%

lambda=[1850,1750,1750,3550]*(1e-9);
Gamma=[ 1.0174e+11, 9.8728e+10,9.7524e+10,1.0174e+11];
mu=[0.187,0.182,0.181,0.349];
gamma_OD=Gamma/(ODN*vol)
simtime = 0:1:1500;



  for  S=1:length(Ein)
      
   [hour, y] = ode45(@(t,y) MRmodel(t,y,mu(S),lambda(S),Ein(S),dep,vol,K,Gamma(S),ODN),simtime, N0);

   OD(:,S)=y./(ODN*vol);
    
  end
     
figure;
plot(hour,OD(:,1),'b-','Linewidth',2);
hold on
plot(hour,OD(:,2),'r-','Linewidth',2);
hold on
plot(hour,OD(:,3),'k-','Linewidth',2);
hold on
plot(hour,OD(:,4),'g-','Linewidth',2);
hold on

plot(time_data,od_data_for_MRmodel(1,:),'bo','Markersize',8);
 hold on
 plot(time_data,od_data_for_MRmodel(2,:),'ro','Markersize',8);
 hold on
 plot(time_data,od_data_for_MRmodel(3,:),'ko','Markersize',8);
 hold on
 plot(time_data,od_data_for_MRmodel(4,:),'go','Markersize',8);
 hold on

     xlabel('Time \it{t} \rm{(hour)}');
    ylabel('OD(\it{t}\rm{)}');
    set(gca, 'FontSize',21);
      xlim([0,1200])
     ylim([0,20])

     
 
%% function that returns dy/dt
 
function dy = MRmodel(t,y,mu,lambda,E0,dep,vol,K,Gamma,ODN)

dy = mu.*(lightpercell(E0,y,dep,vol,K,ODN)/(lambda + lightpercell(E0,y,dep,vol,K,ODN)))*(1-y/Gamma)*y;
        
end


%% function lightpercell
function LPC=lightpercell(L,N,dep,vol,K,ODN)

        cell_conc=@(N) N/vol;
        LPC=L*(1 - 10.^(-K.*cell_conc(N)*dep))./N;

end

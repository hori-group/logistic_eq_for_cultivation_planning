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
 
%% Full model identified parameters
    % Parameter matrix 
    % row: light fulx into culture, column: initial media concentration
    %               C0=1/8 C0=1/4 C0=1/2 C0=1
    %    Ein=0.274
    %    Ein=0.521
    %    Ein=1.09
    %    Ein=2.92
                alpha=[8.7, 8.7, 8.7, 8.7;
                    8.7, 8.7, 8.5, 8.7;
                    8.7, 8.8, 8.7, 8.7;
                    8.7, 8.7, 8.6, 8.7]*1e-12;


                xi=[0.015,0.012,0.013,0.012;
                    0.014,0.013,0.007,0.012;
                    0.013,0.011,0.01,0.014;
                    0.012,0.012,0.011,0.01];
       
%% simulation time   
simtime = 0:1:1500;

%%
%% main code for Ein=0.274
%%

L_ind=1; %at Ein=0.274

for M=1:length(C0)

    y0=[N0,C0(M)];
    [hour, y] = ode23s(@(t,y) Full_model(t,y,alpha(L_ind,M),mu,lambda,Ein(L_ind),dep,vol,K,Gamma,xi(L_ind,M),ODN),...
                                simtime, y0);
   OD(:,M)=y(:,1)./(ODN*vol);
    
end  

   figure
    plot(hour,OD(:,1),'b-','Linewidth',2)
    hold on
    plot(hour,OD(:,2),'r-','Linewidth',2)
    hold on
    plot(hour,OD(:,3),'k-','Linewidth',2)
    hold on
    plot(hour,OD(:,4),'g-','Linewidth',2)
    hold on

    plot(time_data,od_data(1,:),'bo','Markersize',8)
    hold on    
    plot(time_data,od_data(2,:),'ro','Markersize',8)
    hold on
    plot(time_data,od_data(3,:),'ko','Markersize',8)
    hold on
    plot(time_data,od_data(4,:),'go','Markersize',8)
  
    xlabel('Time \it{t} \rm{(hour)}');
    ylabel('OD(\it{t}\rm{)}');
    set(gca, 'FontSize',21);
  
     xlim([0,1200])
     ylim([0,20])

 
%%
%% main code for Ein=0.521
%%

L_ind=2; %at Ein=0.521

for M=1:length(C0)

    y0=[N0,C0(M)];
    [hour, y] = ode23s(@(t,y) Full_model(t,y,alpha(L_ind,M),mu,lambda,Ein(L_ind),dep,vol,K,Gamma,xi(L_ind,M),ODN),...
                                simtime, y0);
   OD(:,M)=y(:,1)./(ODN*vol);
    
end  

   figure
    plot(hour,OD(:,1),'b-','Linewidth',2)
    hold on
    plot(hour,OD(:,2),'r-','Linewidth',2)
    hold on
    plot(hour,OD(:,3),'k-','Linewidth',2)
    hold on
    plot(hour,OD(:,4),'g-','Linewidth',2)
    hold on

    plot(time_data,od_data(5,:),'bo','Markersize',8)
    hold on    
    plot(time_data,od_data(6,:),'ro','Markersize',8)
    hold on
    plot(time_data,od_data(7,:),'ko','Markersize',8)
    hold on
    plot(time_data,od_data(8,:),'go','Markersize',8)
  
    xlabel('Time \it{t} \rm{(hour)}');
    ylabel('OD(\it{t}\rm{)}');
    set(gca, 'FontSize',21);
  
     xlim([0,1200])
     ylim([0,20])
    
%%
%% main code for Ein=1.09
%%

L_ind=3; %at Ein=1.09

for M=1:length(C0)

    y0=[N0,C0(M)];
    [hour, y] = ode23s(@(t,y) Full_model(t,y,alpha(L_ind,M),mu,lambda,Ein(L_ind),dep,vol,K,Gamma,xi(L_ind,M),ODN),...
                                simtime, y0);
   OD(:,M)=y(:,1)./(ODN*vol);
    
end  

   figure
    plot(hour,OD(:,1),'b-','Linewidth',2)
    hold on
    plot(hour,OD(:,2),'r-','Linewidth',2)
    hold on
    plot(hour,OD(:,3),'k-','Linewidth',2)
    hold on
    plot(hour,OD(:,4),'g-','Linewidth',2)
    hold on

    plot(time_data,od_data(9,:),'bo','Markersize',8)
    hold on    
    plot(time_data,od_data(10,:),'ro','Markersize',8)
    hold on
    plot(time_data,od_data(11,:),'ko','Markersize',8)
    hold on
    plot(time_data,od_data(12,:),'go','Markersize',8)
  
    xlabel('Time \it{t} \rm{(hour)}');
    ylabel('OD(\it{t}\rm{)}');
    set(gca, 'FontSize',21);
  
     xlim([0,1200])
     ylim([0,20])

%%
%% main code for Ein=2.92
%%

L_ind=4; %at Ein=2.92

for M=1:length(C0)

    y0=[N0,C0(M)];
    [hour, y] = ode23s(@(t,y) Full_model(t,y,alpha(L_ind,M),mu,lambda,Ein(L_ind),dep,vol,K,Gamma,xi(L_ind,M),ODN),...
                                simtime, y0);
   OD(:,M)=y(:,1)./(ODN*vol);
    
end  

   figure
    plot(hour,OD(:,1),'b-','Linewidth',2)
    hold on
    plot(hour,OD(:,2),'r-','Linewidth',2)
    hold on
    plot(hour,OD(:,3),'k-','Linewidth',2)
    hold on
    plot(hour,OD(:,4),'g-','Linewidth',2)
    hold on

    plot(time_data,od_data(13,:),'bo','Markersize',8)
    hold on    
    plot(time_data,od_data(14,:),'ro','Markersize',8)
    hold on
    plot(time_data,od_data(15,:),'ko','Markersize',8)
    hold on
    plot(time_data,od_data(16,:),'go','Markersize',8)
  
    xlabel('Time \it{t} \rm{(hour)}');
    ylabel('OD(\it{t}\rm{)}');
    set(gca, 'FontSize',21);
  
     xlim([0,1200])
     ylim([0,20])

%% function that returns dy/dt

function dy=Full_model (t,y,alpha,mu,lambda,Ein,dep,vol,K,Gamma,xi,ODN)
 
  dy(1,:)= mu.*(y(2)/(xi+y(2)))*(lightpercell(Ein,y(1),dep,vol,K,ODN)/(lambda + lightpercell(Ein,y(1),dep,vol,K,ODN)))*(1-y(1)/Gamma)*y(1);

  dy(2,:)= -alpha*dy(1,:);
   
  
end


%% function lightpercell
function LPC=lightpercell(L,N,dep,vol,KK,ODN)

        cell_conc=@(N) N/vol;
        LPC=L*(1 - 10.^(-KK.*cell_conc(N)*dep))./N;

end



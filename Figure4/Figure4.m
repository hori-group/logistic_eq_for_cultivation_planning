%%
%% Figure4
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
ein_data=[96.8, 184.4, 386.7, 1034]; %incident light flux density ein (unit is \mu E s^{-1} m^{-2})
light_area=0.002826 %area (unit is m^2)
Ein_data=ein_data*light_area %incident light flux into culture (unit is \mu E s^{-1}) 


time_conversion=3600 %unit is hour s^{-1}

ein=[10:1:100,110:50:1500];

Ein_v=ein*light_area;
Ein=Ein_v*time_conversion


C0=[1/8,1/4,1/2,1]

%% Simulation paramter (see Table2)
 
 mu=0.194;
 lambda=1900*(1e-9);
 Gamma=9.9932e+10;
 Gamma_OD= Gamma/(ODN*vol)
 
alpha=8.7e-12;
xi=0.012;

%% simulation time
simtime = 0:1:25000;

%% 
%% Analytical solution for ODmax
%%

 Nmax_sim=N0+C0/alpha;
 
 ODmax_sim=Nmax_sim/(ODN*vol)
 
 index=find(ODmax_sim>=Gamma_OD);
 
 ODmax_sim(index)=Gamma_OD
 
%% Target OD around 70-90%
upper_rate=0.9
lower_rate=0.7

OD_upp_sim=ODmax_sim*upper_rate
OD_low_sim=ODmax_sim*lower_rate
 

%%
%% Main code
%% C0=0.125

m_ind=1;
od_data_v1=od_data(1:4:13,:);

ind_upp_data=[];
ind_low_data=[];

for i=1:4
 %Experimental data
   if isempty(find(od_data_v1(i,:)>OD_upp_sim(m_ind),1))
       
       time_upp_data_finded(i)=NaN;
       OD_upp_data_finded(i)=NaN;
       time_low_data_finded(i)=NaN;
       OD_low_data_finded(i)=NaN;
       
   else
     %-------------------OD upper---------------------------
     ind_upp_data(i)=find(od_data_v1(i,:)>OD_upp_sim(m_ind),1)
     ind_upp_data(i)=ind_upp_data(i)-1;
   
     OD_upp_data_finded(i)=od_data_v1(i,ind_upp_data(i));
     time_upp_data_finded(i)=time_data(ind_upp_data(i));
    
     %-------------------OD lower-----------------------------
     ind_low_data(i)=find(od_data_v1(i,:)>OD_low_sim(m_ind),1);
     
     OD_low_data_finded(i)=od_data_v1(i,ind_low_data(i));
     time_low_data_finded(i)=time_data(ind_low_data(i));
   end
end


%----Full model-----------------------
%----simulation start-----------------

OD=zeros(length(simtime),length(Ein));

for s=1:length(Ein)
%-------------simulation-----------------------------------------    
    y0=[N0,C0(m_ind)];
    [hour, y] = ode23s(@(t,y) Full_model(t,y,alpha,mu,lambda,Ein(s),dep,vol,K,Gamma,xi,ODN),...
                                simtime, y0);
           
   OD(:,s)=y(:,1)./(ODN*vol);
   
%-------------------OD upper-----------------------------------------------     
     
     ind_upp_sim(s)=find(OD(:,s)>OD_upp_sim(m_ind),1);
     ind_upp_sim(s)=ind_upp_sim(s)-1;   
     
     OD_upp_sim_finded(s)=OD(ind_upp_sim(s),s);
     time_upp_sim_finded(s)=hour(ind_upp_sim(s));

        
     
%-------------------OD lower-----------------------------------------------     
     ind_low_sim(s)=find(OD(:,s)>OD_low_sim(m_ind),1);
     
     OD_low_sim_finded(s)=OD(ind_low_sim(s),s);
     time_low_sim_finded(s)=hour(ind_low_sim(s));
     
     
     
end  

%OD_upp_sim_finded
%OD_low_sim_finded

%time_upp_sim_finded
%time_low_sim_finded

OD_upp_data_finded
time_upp_data_finded

OD_low_data_finded
time_low_data_finded


figure;
p_m_u=plot(Ein_v,time_upp_sim_finded,'b--','Linewidth',1.5);
hold on
p_m_l=plot(Ein_v,time_low_sim_finded,'b-.','Linewidth',1.5);
hold on

p_d_u=plot(Ein_data,time_upp_data_finded,'rv','Markersize',10.5);
hold on
p_d_l=plot(Ein_data,time_low_data_finded,'r^','Markersize',10.5);
hold on

plot([Ein_data(1),Ein_data(1)],[time_low_data_finded(1),time_upp_data_finded(1)],'r-','Linewidth',1.5)
hold on
plot([Ein_data(2),Ein_data(2)],[time_low_data_finded(2),time_upp_data_finded(2)],'r-','Linewidth',1.5)
hold on
plot([Ein_data(3),Ein_data(3)],[time_low_data_finded(3),time_upp_data_finded(3)],'r-','Linewidth',1.5)
hold on
plot([Ein_data(4),Ein_data(4)],[time_low_data_finded(4),time_upp_data_finded(4)],'r-','Linewidth',1.5)
hold on

ar=area(Ein_v',[time_upp_sim_finded' (time_low_sim_finded-time_upp_sim_finded)']); 
set(ar,'Facecolor','None','Linestyle','None','Showbaseline','off');
set(ar(1),'Facecolor','None','Linestyle','None','Showbaseline','off');
set(ar(2),'Facecolor',[0,0.2,1.0],'FaceAlpha',0.2,'Linestyle','None');

ylim([0,1200])
xlim([0,3.5])
   set(gca, 'FontSize',23);    


%%
%% Main code
%% C0=0.25

m_ind=2;
od_data_v1=od_data(2:4:14,:);

ind_upp_data=[];
ind_low_data=[];

for i=1:4
 %experimnetal data
   if isempty(find(od_data_v1(i,:)>OD_upp_sim(m_ind),1))
       
       time_upp_data_finded(i)=NaN;
       OD_upp_data_finded(i)=NaN;
       time_low_data_finded(i)=NaN;
       OD_low_data_finded(i)=NaN;
       
   else
     %-------------------OD upper---------------------------
     ind_upp_data(i)=find(od_data_v1(i,:)>OD_upp_sim(m_ind),1)
     ind_upp_data(i)=ind_upp_data(i)-1;
   
     OD_upp_data_finded(i)=od_data_v1(i,ind_upp_data(i));
     time_upp_data_finded(i)=time_data(ind_upp_data(i));
    
     %-------------------OD lower-----------------------------
     ind_low_data(i)=find(od_data_v1(i,:)>OD_low_sim(m_ind),1);
     
     OD_low_data_finded(i)=od_data_v1(i,ind_low_data(i));
     time_low_data_finded(i)=time_data(ind_low_data(i));
   end
end


%----Full model-----------------------
%----simulation start-----------------

OD=zeros(length(simtime),length(Ein));

for s=1:length(Ein)
%-------------simulation-----------------------------------------    
    y0=[N0,C0(m_ind)];
    [hour, y] = ode23s(@(t,y) Full_model(t,y,alpha,mu,lambda,Ein(s),dep,vol,K,Gamma,xi,ODN),...
                                simtime, y0);
           
   OD(:,s)=y(:,1)./(ODN*vol);
   
%-------------------OD upper-----------------------------------------------     
     
     ind_upp_sim(s)=find(OD(:,s)>OD_upp_sim(m_ind),1);
     ind_upp_sim(s)=ind_upp_sim(s)-1;   
     
     OD_upp_sim_finded(s)=OD(ind_upp_sim(s),s);
     time_upp_sim_finded(s)=hour(ind_upp_sim(s));

        
     
%-------------------OD lower-----------------------------------------------     
     ind_low_sim(s)=find(OD(:,s)>OD_low_sim(m_ind),1);
     
     OD_low_sim_finded(s)=OD(ind_low_sim(s),s);
     time_low_sim_finded(s)=hour(ind_low_sim(s));
     
     
     
end  

%OD_upp_sim_finded
%OD_low_sim_finded

%time_upp_sim_finded
%time_low_sim_finded

OD_upp_data_finded
time_upp_data_finded

OD_low_data_finded
time_low_data_finded


figure;
p_m_u=plot(Ein_v,time_upp_sim_finded,'b--','Linewidth',1.5);
hold on
p_m_l=plot(Ein_v,time_low_sim_finded,'b-.','Linewidth',1.5);
hold on

p_d_u=plot(Ein_data,time_upp_data_finded,'rv','Markersize',10.5);
hold on
p_d_l=plot(Ein_data,time_low_data_finded,'r^','Markersize',10.5);
hold on

plot([Ein_data(1),Ein_data(1)],[time_low_data_finded(1),time_upp_data_finded(1)],'r-','Linewidth',1.5)
hold on
plot([Ein_data(2),Ein_data(2)],[time_low_data_finded(2),time_upp_data_finded(2)],'r-','Linewidth',1.5)
hold on
plot([Ein_data(3),Ein_data(3)],[time_low_data_finded(3),time_upp_data_finded(3)],'r-','Linewidth',1.5)
hold on
plot([Ein_data(4),Ein_data(4)],[time_low_data_finded(4),time_upp_data_finded(4)],'r-','Linewidth',1.5)
hold on

ar=area(Ein_v',[time_upp_sim_finded' (time_low_sim_finded-time_upp_sim_finded)']); 
set(ar,'Facecolor','None','Linestyle','None','Showbaseline','off');
set(ar(1),'Facecolor','None','Linestyle','None','Showbaseline','off');
set(ar(2),'Facecolor',[0,0.2,1.0],'FaceAlpha',0.2,'Linestyle','None');

ylim([0,1200])
xlim([0,3.5])
   set(gca, 'FontSize',23);    

   
%%
%% Main code
%% C0=0.5

m_ind=3;
od_data_v1=od_data(3:4:15,:);

ind_upp_data=[];
ind_low_data=[];

for i=1:4
 %experimnetal data
   if isempty(find(od_data_v1(i,:)>OD_upp_sim(m_ind),1))
       
       time_upp_data_finded(i)=NaN;
       OD_upp_data_finded(i)=NaN;
       time_low_data_finded(i)=NaN;
       OD_low_data_finded(i)=NaN;
       
   else
     %-------------------OD upper---------------------------
     ind_upp_data(i)=find(od_data_v1(i,:)>OD_upp_sim(m_ind),1)
     ind_upp_data(i)=ind_upp_data(i)-1;
   
     OD_upp_data_finded(i)=od_data_v1(i,ind_upp_data(i));
     time_upp_data_finded(i)=time_data(ind_upp_data(i));
    
     %-------------------OD lower-----------------------------
     ind_low_data(i)=find(od_data_v1(i,:)>OD_low_sim(m_ind),1);
     
     OD_low_data_finded(i)=od_data_v1(i,ind_low_data(i));
     time_low_data_finded(i)=time_data(ind_low_data(i));
   end
end


%----Full model-----------------------
%----simulation start-----------------

OD=zeros(length(simtime),length(Ein));

for s=1:length(Ein)
%-------------simulation-----------------------------------------    
    y0=[N0,C0(m_ind)];
    [hour, y] = ode23s(@(t,y) Full_model(t,y,alpha,mu,lambda,Ein(s),dep,vol,K,Gamma,xi,ODN),...
                                simtime, y0);
           
   OD(:,s)=y(:,1)./(ODN*vol);
   
%-------------------OD upper-----------------------------------------------     
     
     ind_upp_sim(s)=find(OD(:,s)>OD_upp_sim(m_ind),1);
     ind_upp_sim(s)=ind_upp_sim(s)-1;   
     
     OD_upp_sim_finded(s)=OD(ind_upp_sim(s),s);
     time_upp_sim_finded(s)=hour(ind_upp_sim(s));

        
     
%-------------------OD lower-----------------------------------------------     
     ind_low_sim(s)=find(OD(:,s)>OD_low_sim(m_ind),1);
     
     OD_low_sim_finded(s)=OD(ind_low_sim(s),s);
     time_low_sim_finded(s)=hour(ind_low_sim(s));
     
     
     
end  

%OD_upp_sim_finded
%OD_low_sim_finded

%time_upp_sim_finded
%time_low_sim_finded

OD_upp_data_finded
time_upp_data_finded

OD_low_data_finded
time_low_data_finded


figure;
p_m_u=plot(Ein_v,time_upp_sim_finded,'b--','Linewidth',1.5);
hold on
p_m_l=plot(Ein_v,time_low_sim_finded,'b-.','Linewidth',1.5);
hold on

p_d_u=plot(Ein_data,time_upp_data_finded,'rv','Markersize',10.5);
hold on
p_d_l=plot(Ein_data,time_low_data_finded,'r^','Markersize',10.5);
hold on

plot([Ein_data(1),Ein_data(1)],[time_low_data_finded(1),time_upp_data_finded(1)],'r-','Linewidth',1.5)
hold on
plot([Ein_data(2),Ein_data(2)],[time_low_data_finded(2),time_upp_data_finded(2)],'r-','Linewidth',1.5)
hold on
plot([Ein_data(3),Ein_data(3)],[time_low_data_finded(3),time_upp_data_finded(3)],'r-','Linewidth',1.5)
hold on
plot([Ein_data(4),Ein_data(4)],[time_low_data_finded(4),time_upp_data_finded(4)],'r-','Linewidth',1.5)
hold on

ar=area(Ein_v',[time_upp_sim_finded' (time_low_sim_finded-time_upp_sim_finded)']); 
set(ar,'Facecolor','None','Linestyle','None','Showbaseline','off');
set(ar(1),'Facecolor','None','Linestyle','None','Showbaseline','off');
set(ar(2),'Facecolor',[0,0.2,1.0],'FaceAlpha',0.2,'Linestyle','None');

ylim([0,1200])
xlim([0,3.5])
   set(gca, 'FontSize',23);    


%%
%% Main code
%% C0=1
m_ind=4;
od_data_v1=od_data(4:4:16,:);

ind_upp_data=[];
ind_low_data=[];

for i=1:4
 %experimnetal data
   if isempty(find(od_data_v1(i,:)>OD_upp_sim(m_ind),1))
       
       time_upp_data_finded(i)=NaN;
       OD_upp_data_finded(i)=NaN;
       time_low_data_finded(i)=NaN;
       OD_low_data_finded(i)=NaN;
       
   else
     %-------------------OD upper---------------------------
     ind_upp_data(i)=find(od_data_v1(i,:)>OD_upp_sim(m_ind),1)
     ind_upp_data(i)=ind_upp_data(i)-1;
   
     OD_upp_data_finded(i)=od_data_v1(i,ind_upp_data(i));
     time_upp_data_finded(i)=time_data(ind_upp_data(i));
    
     %-------------------OD lower-----------------------------
     ind_low_data(i)=find(od_data_v1(i,:)>OD_low_sim(m_ind),1);
     
     OD_low_data_finded(i)=od_data_v1(i,ind_low_data(i));
     time_low_data_finded(i)=time_data(ind_low_data(i));
   end
end


%----Full model-----------------------
%----simulation start-----------------

OD=zeros(length(simtime),length(Ein));

for s=1:length(Ein)
%-------------simulation-----------------------------------------    
    y0=[N0,C0(m_ind)];
    [hour, y] = ode23s(@(t,y) Full_model(t,y,alpha,mu,lambda,Ein(s),dep,vol,K,Gamma,xi,ODN),...
                                simtime, y0);
           
   OD(:,s)=y(:,1)./(ODN*vol);
   
%-------------------OD upper-----------------------------------------------     
     
     ind_upp_sim(s)=find(OD(:,s)>OD_upp_sim(m_ind),1);
     ind_upp_sim(s)=ind_upp_sim(s)-1;   
     
     OD_upp_sim_finded(s)=OD(ind_upp_sim(s),s);
     time_upp_sim_finded(s)=hour(ind_upp_sim(s));

        
     
%-------------------OD lower-----------------------------------------------     
     ind_low_sim(s)=find(OD(:,s)>OD_low_sim(m_ind),1);
     
     OD_low_sim_finded(s)=OD(ind_low_sim(s),s);
     time_low_sim_finded(s)=hour(ind_low_sim(s));
     
     
     
end  

%OD_upp_sim_finded
%OD_low_sim_finded

%time_upp_sim_finded
%time_low_sim_finded

OD_upp_data_finded
time_upp_data_finded

OD_low_data_finded
time_low_data_finded


figure;
p_m_u=plot(Ein_v,time_upp_sim_finded,'b--','Linewidth',1.5);
hold on
p_m_l=plot(Ein_v,time_low_sim_finded,'b-.','Linewidth',1.5);
hold on

p_d_u=plot(Ein_data,time_upp_data_finded,'rv','Markersize',10.5);
hold on
p_d_l=plot(Ein_data,time_low_data_finded,'r^','Markersize',10.5);
hold on

plot([Ein_data(1),Ein_data(1)],[time_low_data_finded(1),time_upp_data_finded(1)],'r-','Linewidth',1.5)
hold on
plot([Ein_data(2),Ein_data(2)],[time_low_data_finded(2),time_upp_data_finded(2)],'r-','Linewidth',1.5)
hold on
plot([Ein_data(3),Ein_data(3)],[time_low_data_finded(3),time_upp_data_finded(3)],'r-','Linewidth',1.5)
hold on
plot([Ein_data(4),Ein_data(4)],[time_low_data_finded(4),time_upp_data_finded(4)],'r-','Linewidth',1.5)
hold on

ar=area(Ein_v',[time_upp_sim_finded' (time_low_sim_finded-time_upp_sim_finded)']); 
set(ar,'Facecolor','None','Linestyle','None','Showbaseline','off');
set(ar(1),'Facecolor','None','Linestyle','None','Showbaseline','off');
set(ar(2),'Facecolor',[0,0.2,1.0],'FaceAlpha',0.2,'Linestyle','None');

ylim([0,1200])
xlim([0,3.5])
   set(gca, 'FontSize',23);    

     

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




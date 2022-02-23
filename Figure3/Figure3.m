%%
%% Figure3
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

    ODmax_data=[max(od_data(1,:)),max(od_data(2,:)),max(od_data(3,:)),max(od_data(4,:));
        max(od_data(5,:)),max(od_data(6,:)),max(od_data(7,:)),max(od_data(8,:));
        max(od_data(9,:)),max(od_data(10,:)),max(od_data(11,:)),max(od_data(12,:));
        max(od_data(13,:)),max(od_data(14,:)),max(od_data(15,:)),max(od_data(16,:))]
    
    %Three experimental data is not reached before stationary phase.
    %These cannot be searched for ODmax
    ODmax_data(1,3)=NaN %Ein=0.274 C0=0.5
    ODmax_data(1,4)=NaN %Ein=0.274 C0=1
    ODmax_data(2,4)=NaN %Ein=0.521 C0=1

%% initial values
 ODN=30.10*1e6;
 vol = 200; %Volume V (unit is mL)
 
 OD0=0.025 %initial OD
 N0 = OD0*ODN*vol; %initial OD fixed

%% initial media concentration

C0=[0:0.01:1.5];

%% Identified paramter (see Table2) for ODmax
 
 Gamma=9.9932e+10;
 Gamma_OD= Gamma/(ODN*vol)

 alpha=8.7e-12;

%% 
%% Main code for analytical solution for ODmax
%%

 Nmax_sim_C_depend=N0+C0/alpha;
 
 ODmax_sim_C_depend=Nmax_sim_C_depend/(ODN*vol)
 
 ODmax_sim=ODmax_sim_C_depend;
 index=find(ODmax_sim>=Gamma_OD);
 
  ODmax_sim(index)=Gamma_OD
  
  
  %% figure
  figure
  m=plot(C0,ODmax_sim,'k-','Linewidth',2);
  hold on
  l1=plot([0.125,0.25,0.5,1],ODmax_data(1,:),'bo','Markersize',10);
  hold on
  l2= plot([0.125,0.25,0.5,1],ODmax_data(2,:),'rs','Markersize',10);
   hold on
   l3=plot([0.125,0.25,0.5,1],ODmax_data(3,:),'md','Markersize',10);
   hold on
   l4= plot([0.125,0.25,0.5,1],ODmax_data(4,:),'gh','Markersize',10);
         
   ylabel({'Maximum OD'});
   xlabel('Initial media concentration');
      ylim([0,20])
      xlim([0,1.5])
      set(gca, 'FontSize',21); 
    
  


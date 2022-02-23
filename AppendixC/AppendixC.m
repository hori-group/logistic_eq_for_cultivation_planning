%%
%% AppendixC
%%

clear; 
close all;

%% excel file
    filename = 'AppendixC.xlsx';
 
    sheet=1; 
    xlRange='b18:j42';
    d_in=xlsread(filename,sheet,xlRange);
    
 Depth=d_in(2:8,1);
 
 OD=d_in(1,2:7);
 
 log_eout_per_ein=d_in(2:length(Depth)+1,2:length(OD)+1);
 
 K_OD=d_in(14:19,2:3) %Row1: OD Row2: KOD
 K_D=d_in(20:25,2:3)  %Row1: D Row2: KD
 
 
 %% Panel A
 figure;
 plot(Depth,log_eout_per_ein(:,1),'o','Markersize',8);
 hold on
 plot(Depth,log_eout_per_ein(:,2),'o','Markersize',8);
 hold on
 plot(Depth,log_eout_per_ein(:,3),'o','Markersize',8);
 hold on
 plot(Depth,log_eout_per_ein(:,4),'o','Markersize',8);
 hold on
 plot(Depth,log_eout_per_ein(:,5),'o','Markersize',8);
 hold on
 plot(Depth,log_eout_per_ein(:,6),'o','Markersize',8);
hold on


 plot(Depth,K_D(1,2)*Depth,'-','Linewidth',2);
 hold on
 plot(Depth,K_D(2,2)*Depth,'-','Linewidth',2);
 hold on
plot(Depth,K_D(3,2)*Depth,'-','Linewidth',2);
 hold on
plot(Depth,K_D(4,2)*Depth,'-','Linewidth',2);
 hold on
 plot(Depth,K_D(5,2)*Depth,'-','Linewidth',2);
 hold on
 plot(Depth,K_D(6,2)*Depth,'-','Linewidth',2);
 hold on
 
   newcolors=[0.64,0.08,0.18; 
     0.00,0.45,0.74;
     0.85,0.33,0.10;
     0.93,0.69,0.13;
     0.49,0.18,0.56;
     0.47,0.67,0.19];
 colororder(newcolors);
box on
xlabel('Depth \it{D} \rm{(cm)}');
    xlim([0,6])
    ylim([-3,0])
    lgd=legend('0.149','0.31','0.6','1.17','2.47','4.78');
    lgd.NumColumns=2;
    title(lgd,['OD'])
    lgd.FontSize = 18;
    lgd.Location='Southwest';
    set(gca, 'FontSize',23);
    
    
    
%% Panel B
 figure;
 plot(OD,log_eout_per_ein(2,:),'o','Markersize',8);
 hold on
 plot(OD,log_eout_per_ein(3,:),'o','Markersize',8);
 hold on
 plot(OD,log_eout_per_ein(4,:),'o','Markersize',8);
 hold on
 plot(OD,log_eout_per_ein(5,:),'o','Markersize',8);
 hold on
 plot(OD,log_eout_per_ein(6,:),'o','Markersize',8);
hold on
 plot(OD,log_eout_per_ein(7,:),'o','Markersize',8);
hold on

OD_sim_value=[0,OD,5,6];
 plot(OD_sim_value,K_OD(1,2)*OD_sim_value,'-','Linewidth',2);
 hold on
plot(OD_sim_value,K_OD(2,2)*OD_sim_value,'-','Linewidth',2);
 hold on
plot(OD_sim_value,K_OD(3,2)*OD_sim_value,'-','Linewidth',2);
 hold on
 plot(OD_sim_value,K_OD(4,2)*OD_sim_value,'-','Linewidth',2);
 hold on
 plot(OD_sim_value,K_OD(5,2)*OD_sim_value,'-','Linewidth',2);
 hold on
 plot(OD_sim_value,K_OD(6,2)*OD_sim_value,'-','Linewidth',2);
 hold on
 
newcolors=[0.64,0.08,0.18; 
     0.00,0.45,0.74;
     0.85,0.33,0.10;
     0.93,0.69,0.13;
     0.49,0.18,0.56;
     0.47,0.67,0.19];
  colororder(newcolors);
box on
xlabel('OD');
    xlim([0,6])
    ylim([-3,0])
    lgd=legend('1','2','3','4','5','6');
    lgd.NumColumns=2;
    title(lgd,['\it{D}\rm{ (cm)}'])
    lgd.FontSize = 18;
    lgd.Location='Southwest';
    set(gca, 'FontSize',23);
    
      
     
%% Panel C
figure;
plot(K_OD(:,1),K_OD(:,2),'bs','Markersize',10)
hold on
plot(K_D(:,1),K_D(:,2),'bd','Markersize',10)
hold on
x=[0:6];
plot(x,-0.1529*x,'b-','Linewidth',2);
 xlabel('Depth \it{D}\rm{, OD}');
    lgd.FontSize = 23;
    lgd.Location='Southwest';
    set(gca, 'FontSize',23);
    xlim([0,6])
    ylim([-1,0])

 
clear all
format long 
clear all
close all
clc




%% 
T = readtable('household_power_consumption.txt',"Delimiter",';');
%fid=fopen('household_power_consumption.txt');
%f=textscan(fid,'%s', 'delimiter','\n');
%fclose(fid);
%figure
    
%%
% %reading the data for first and third rows 
%D is the variable for active power consumption

D =table2array(T(1:11949,2));
%E is the Global reactive power
E=table2array(T(1:11949,3))
%F is voltage 
F= table2array(T(1:11949,4)); 

%G is the global intensity 
G = table2array(T(1:11949,5)); 


date=D';
% %making transpose
S =table2array(T(1:11949,3));
global_active_power = S';


global_reactive_power=table2array(T(1:11949,4));
global_reactive_power_1 = global_reactive_power';

voltage = table2array(T(1:11949,5));
voltage_1= voltage';

global_intensity = table2array(T(1:11949,6));
global_intensity_1= global_intensity';




%%
%making a 2 row matrix from it
%f = [date; global_e]
% 
% 
%%
%plotting the data for various varaibles
figure
plot(date,global_active_power,"r-");
hold on;
plot(date,global_reactive_power_1,"b-");
xlabel('Date','FontSize',14)
ylabel('Total Energy Consumption/Total reactive power ','FontSize',14);
legend('Global Reactive Power','Total Reactive Power')
%xlabel('date','Font',1

figure
plot(date,voltage_1,"g-");
xlabel('Date','FontSize',14)
ylabel('Voltage w.r.t. to date','FontSize',14);




figure 
plot(date, global_intensity_1, "c-")
xlabel('Date','FontSize',14)
ylabel('Global Intensity','FontSize',14);


%
% %%
% %making data at an interva; of 
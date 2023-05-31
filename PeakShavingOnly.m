clear
close all

clc

load('Price')
load('load_data.mat')

hour = 1;
ts = 4/3600;
T = hour*1/ts;
tt = 4:4:3600*hour;

P_g = avg_active_grid(1:T,1);

e_ch=0.95;
e_dis =0.95;
del_max = 2000;
del_min = -del_max;
B_0 = 1000*ones(T,1);
B_max = 2000*ones(T,1);
B_min = 200*ones(T,1);
h=0.25;
S_B_max = 1.5*del_max*h/e_ch;
x_upper= del_max*h*ones(T,1);
x_lower= del_min*h*ones(T,1);
q_upper= S_B_max*ones(T,1);
q_lower= -S_B_max*ones(T,1);

C_cap =S_B_max^2*ones(T,1);
M = tril(ones(T,T));
ub=1e6;
lb=1;
%%
battery.n = 5000;
battery.cell = 0.3;% Cell price, $/Wh
battery.energy = del_max*(3/60); % MWh, battery energy, up for 3min
battery.socmax = 0.8; 
battery.socmin = 0.2;
battery.socini = 0.6;

lambda.battery = battery.cell*10^(6)/...
    (2*battery.n*(battery.socmax-battery.socmin)); %Battery cost, $/MWh
 




%% Convex optimization formulation
%E=sum(A,2);
%P_g = sum(A,2).*ones(T,1);
%Q_g = sum(A,2)*0.3.*ones(T,1);
C_cap =S_B_max^2.*ones(T,1);
M = tril(ones(T,T));
ub=1e6;
lb=1;
ts=0.0011;
cvx_begin 

variables x_ch(T,1) x_ds(T,1)  
variables b_p(T,1) b_p_bat(T,1)
minimize sum (100*Price.*(P_g-b_p)*ts)+ lambda.battery*norm(b_p,1)*ts% 
subject to 
    zeros(T,1)<= x_ch <= x_upper;    %%charging 
    zeros(T,1)<= x_ds <= -x_lower;    %% discharging 
    b_p_bat == x_ch-x_ds;               %%change in charge level of battery
    b_p == x_ch/e_ch-x_ds*e_dis;        %%active power output of battery
    B_min <= B_0 + M*b_p_bat <= B_max;  %%battery capacity
   
cvx_end

a=(P_g)-b_p;

B =  B_0+ b_p_bat*h;
%%
    SoC_3 = zeros(T,1);
    SoC_3(1) = battery.socini;
    for i = 2:T
        SoC_3(i) = (battery.socini*battery.energy-sum(b_p(1:i-1))*h)/battery.energy;
    end


%% Bill
total =sum(100*Price.*a*ts) + lambda.battery*norm(b_p,1)*ts

%% plots
figure;

plot(tt,(P_g),tt,a)
legend('Original active power ','Peak shaving only ')
xlabel('Time (s)')
ylabel('Power(W)')

figure;
subplot(3,1,1)
plot(tt,(P_g),tt,a)
legend('Original active power ','Peak shaving only ')
xlabel('Time [s]')
ylabel('Power(W)')
 subplot(3,1,2)
 plot(tt,B,'b');
xlabel('Time[s]');
ylabel('State of battery');
 subplot(3,1,3)
 plot(tt,Price,'b');
xlabel('Time[s]');
ylabel('$');
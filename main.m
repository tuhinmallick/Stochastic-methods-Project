%%==========================================================
% Analysis of eletircity consumption for data for Paris 
% and comparison of regularization methods with GNN 
%%==========================================================
clear all
close all
clc
rand('seed',1)
randn('seed',1)

%%
% Load the data
load energy.mat
all_var=[global_active_power;global_intensity';global_reactive_power']';

%%
%global_Active_power
temp_1= all_var(:,1)';
%global_Reactive_power
temp_2 = all_var(:,2)';
%global_Intensity 
temp_3=all_var(:,3)';

T=length(temp_3);
temp=temp_3;
%% Plot 

%bar(all_var(:,1))
figure

subplot(3,1,1)
plot(temp_1)
ylabel('Global Reactive Power','FontSize',14);

subplot(3,1,2)
ylabel('Global Active Power','FontSize',14);

plot(temp)
ylabel('Voltage ','FontSize',14);


subplot(3,1,3)
plot(temp)
ylabel('Global Reactive power ','FontSize',14);

%% 
% Compute an expectation for every day of the year

expectation_temp=zeros(1,365);
for ind=1:365
    expectation_temp(ind)=sum(temp(ind:365:T))./length(temp(ind:365:T));
    date(ind)=ind;
end
%%
% Form the seasonal trend vector out of expectations for every day
SeasonalTrend=zeros(1,T);
for t=1:T
    ind=mod(t,365);
    if ind==0
        ind=365;
    end
    SeasonalTrend(t)=expectation_temp(ind);
    temp(t)=temp(t)-SeasonalTrend(t);
end

%%
% Create higher order polynomial trends (t^2, t^3, ...,t^P)
% and normalize it (by scaling all of the trends
% to the interval [0 1])
P=40;
figure;plot(temp)
ylabel('Global Voltage','FontSize',14);
%Names={'Constant', 'Seasonal'};
X_trend=zeros(T,P);
for p=0:P
   X(:,p+1)=(1/(T^p)).*[1:T]'.^p; 
   Names{p+1}=['t^{' num2str(p) '}'];
end

 Y=temp';
 
%%
% Solving the original (not regularized) regression problem

min(eig(cov(X)));
Beta_3=ComputeBeta(X,Y);
figure;plot(Beta_3,'k-o','LineWidth',2,'MarkerSize',12);set(gca,'XTick',1:length(Names),'XTickLabel',Names)
title('Global Voltage ')
ylabel('Parameter Values','FontSize',16);
set(gca,'XScale','log','LineWidth',2,'FontSize',16)
axis tight
Y_ML=Beta_3'*X';
%%
% Demonstrate the problem with predictions for D days
D=3000;
X_pred=zeros(T+D,P+1);
for p=0:P
   X_pred(:,p+1)=(1/(T^(p))).*[1:T+D]'.^(p); 
end
title('Global Reactive power')
ff=figure;plot(Y,':.','LineWidth',0.5);hold on;
plot(Beta_3'*X_pred','m.-','LineWidth',2);


%%
% Solve regularized problem for different alpha

alpha=[1e-5 1e-4 1e-3 5e-3 1e-2 5e-2 1e-1 5e-1 1 2 4 8 16 32 64 128 ];
figure;hold on;
tic;
for i=1:length(alpha)
    BetaReg_3(:,i)=ComputeBetaRegularized(X,Y,alpha(i));
    Y_MLReg(i,:)=BetaReg_3(:,i)'*X';
    ErrorNorm(i)=(Y'-Y_MLReg(i,:))*(Y'-Y_MLReg(i,:))';
    SolutionNorm(i)=BetaReg_3(:,i)'*BetaReg_3(:,i);
    plot(SolutionNorm(i),ErrorNorm(i),'o','LineWidth',2);
    title('Beta values')
    text(SolutionNorm(i)+1,ErrorNorm(i)-1,['\alpha=' num2str(alpha(i))],'FontSize',12);
end
time_l2=toc;
plot(SolutionNorm,ErrorNorm,':','LineWidth',2);
xlabel('Norm of the Solution \beta^\alpha','FontSize',16);
ylabel('Norm of the Residuals \epsilon','FontSize',16);
set(gca,'XScale','log','LineWidth',2,'FontSize',16)
axis tight

%%
% Visualize Parameters and identified trend
i=11;
gg=figure;plot(0*BetaReg_3(:,i),'k:','LineWidth',2);hold on;plot(BetaReg_3(:,i),'r-o','LineWidth',2,'MarkerSize',12);set(gca,'XTick',1:length(Names),'XTickLabel',Names)


%%
% Demonstrate the predictions for D days
X_pred=zeros(T+D,P+1);
for p=0:P
   X_pred(:,p+1)=(1/(T^(p))).*[1:T+D]'.^(p); 
end
figure(ff);
plot(BetaReg_3(:,i)'*X_pred','r.-','LineWidth',2);

%%
% Solve l1-regularized problem (lasso regression) 
tic;
[B,FitInfo] = lasso(X(:,2:size(X,2)),Y,'CV',10);
time_l1=toc;

%%
disp(['Time for L2-regularization is ' num2str(time_l2) '; time for l1-regularization is ' num2str(time_l1) '; L2 is ' num2str(time_l1/time_l2) ' times cheaper than L1'])
figure(ff); 
plot(FitInfo.Intercept(FitInfo.IndexMinMSE)+B(1:size(B,1),FitInfo.IndexMinMSE)'*X_pred(:,[2:size(X,2)])','g.-','LineWidth',2);
set(gca,'LineWidth',2,'FontSize',16);
ylabel('Temperature deviations from the seasonal mean','FontSize',14);
ylim([-10 15])
legend('historical data','no regularization','L2-regularized Fit','L1-regularized Fit')
figure(gg);plot([FitInfo.Intercept(FitInfo.IndexMinMSE) B(:,FitInfo.IndexMinMSE)'],'g-o','LineWidth',2,'MarkerSize',9);set(gca,'XTick',1:length(Names),'XTickLabel',Names)
set(gca,'LineWidth',2,'FontSize',16);set(gcf,'Position',[10 100 800  600])
legend({'zero line','L2-regularization','L1-regularization'})
ylabel('Parameter Values','FontSize',16);
set(gca,'XTick',1:length(Names),'XTickLabel',Names)

%%
% We take all the beta  values into a single beta value  by a simple mean
% operation 
load("beta_1.mat")
load("beta_2.mat")
load("beta_3.mat")

total_beta=[Beta_1 Beta_2 Beta_3];

for nn = 1:length(total_beta) 
    final_beta(nn,1) = mean(total_beta(nn,:))
end

%%
% Now using this beta values for prediction 
% We do the prediction of the system using final_beta values which gives

D=3000;
X_pred=zeros(T+D,P+1);
for p=0:P
   X_pred(:,p+1)=(1/(T^(p))).*[1:T+D]'.^(p); 
end
title('Global Reactive power')
ff=figure;plot(Y,':.','LineWidth',0.5);hold on;
plot(final_beta'*X_pred','m.-','LineWidth',2);

%%
C= cov(Beta) 

%%
%making subplots using images
% Load saved figures
c1=hgload('1_p10.fig');
c2=hgload('2_p10.fig');
c3=hgload('3_p10.fig');
k1=hgload('1_p40.fig');
k2=hgload('2_p40.fig');
k3=hgload('3_p40.fig');

% Prepare subplots
figure
h(1)=subplot(2,3,1);
h(2)=subplot(2,3,2);
h(3)=subplot(2,3,3);
h(4)=subplot(2,3,4);
h(5)=subplot(2,3,5);
h(6)=subplot(2,3,6);

% Paste figures on the subplots
copyobj(allchild(get(c1,'CurrentAxes')),h(1));
copyobj(allchild(get(c2,'CurrentAxes')),h(2));
copyobj(allchild(get(c3,'CurrentAxes')),h(3));
copyobj(allchild(get(k1,'CurrentAxes')),h(4));
copyobj(allchild(get(k2,'CurrentAxes')),h(5));
copyobj(allchild(get(k3,'CurrentAxes')),h(6));

% Add legends
l(1)=legend(h(1),'Beta Values P =10')
l(2)=legend(h(2),'Beta Values P =40')

%%%making subplots using images
% Load saved figures
c1=hgload('temp_1.fig');
c2=hgload('temp_2.fig');
c3=hgload('temp_3.fig');

% Prepare subplots
figure
h(1)=subplot(1,3,1);
h(2)=subplot(1,3,2);
h(3)=subplot(1,3,3);

% Paste figures on the subplots
copyobj(allchild(get(c1,'CurrentAxes')),h(1));
copyobj(allchild(get(c2,'CurrentAxes')),h(2));
copyobj(allchild(get(c3,'CurrentAxes')),h(3));
% Add legends
l(1)=legend(h(1),'Beta Values P =10')
l(2)=legend(h(2),'Beta Values P =40')

%%
rse= 
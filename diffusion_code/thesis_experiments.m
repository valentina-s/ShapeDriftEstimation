% This script stores the experiments for simulating diffusions and estimating their parameters

% setting line smoothing on - is it possible in matlab R2012a?
% set(0, 'DefaultLineSmoothing', 'on');


s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s);


% save plots
saveSimulations = true;
saveEstimations = true;

results_path = '../diffusion_results/';

% Diffusion Simulation

% ------------------------- create shapes ---------------------------------

% circle
a = [0:0.1:2*pi];
x = 10*[cos(a)' sin(a)'] + 30;
x_small = x(1:7:end,:);

% ellipse
ellipse_lr = [10*sin(a)' 5*cos(a)'] + 30;
ellipse_ud = [5*sin(a)' 10*cos(a)'] + 30;

% dumbbell
a = 0.9;
y = [-1:0.05:1];
yy_plus = (y.^4.*(a - y.^2)/(a^2));
yy_minus = -y.^4.*(a - y.^2)/(a^2);
dumbbell = [[y';y(end:-1:1)'] [yy_plus'+0.15 ; yy_minus'-0.15]]*15+30;

% ------------------------------- no drift --------------------------------

s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s);


dt = 0.05;
T  = 30;

[bdryPts_t, ctrlPts_t] = Diffusion(x,x_small,zeros(size(x_small)),10,dt,T);
plot3D(bdryPts_t,ctrlPts_t,dt)
title('Brownian Motion','FontSize',14, 'fontweight','bold')

if (saveSimulations == true)
    saveas(gcf,fullfile(results_path,'BrownianMotion1.fig'))
    saveas(gcf,fullfile(results_path,'BrownianMotion1.png'))
end

pause(5)

view(3)

if (saveSimulations == true)
    saveas(gcf,fullfile(results_path,'BrownianMotion2.fig'))
    saveas(gcf,fullfile(results_path,'BrownianMotion2.png'))
end
pause(5)



% ----------------------------- constant drift ----------------------------
% Notes: make it more drifty

s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s);


dt = 0.001;
T = 30;

% theta_constant = 0.1*randn(size(x_small));
theta_constant = 0.5*ones(size(x_small));
[bdryPts_t, ctrlPts_t] = Diffusion(x,x_small,theta_constant,10,dt,T);

plot3D(bdryPts_t,ctrlPts_t,dt)
title('Constant Drift','FontSize',14,'fontweight','bold')

if (saveSimulations == true)
    saveas(gcf,fullfile(results_path,'ConstantDrift1.fig'))
    saveas(gcf,fullfile(results_path,'ConstantDrift1.png'))
end
pause(5)

view(3)
if (saveSimulations == true)
    saveas(gcf,fullfile(results_path,'ConstantDrift2.fig'))
    saveas(gcf,fullfile(results_path,'ConstantDrift2.png'))
end
pause(5)

% -------------------------------- OU drift -------------------------------

s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s);


dt = 0.05; % not sure of 0.05 or 0.01 is better
T = 30;

theta_OU = 0.5;
theta_OU = 1;
%[bdryPts_t ctrlPts_t alpha_t] = Diffusion_drift_OU(1.5*ellipse_lr,1.5*ellipse_lr(1:5:end-4,:),1.5*dumbbell,zeros(100,100),10,10,theta_OU,dt,T);
%[bdryPts_t ctrlPts_t alpha_t] = Diffusion_drift_OU(x,x_small,1.5*dumbbell,zeros(100,100),10,10,theta_OU,dt,T);

[bdryPts_t ctrlPts_t alpha_t] = Diffusion_drift_OU(x+15,x_small+15,1.5*dumbbell,zeros(100,100),10,10,1,dt,T);

plot3D(bdryPts_t,ctrlPts_t,dt)
title('Ornstein-Uhlenbeck Drift','FontSize',14,'fontweight','bold')

% plotting the dumbbell (adding the first point to the end to connect the contour)
plot3((T+1)*ones(size(dumbbell,1)+1,1),1.5*[dumbbell(:,1);dumbbell(1,1)],1.5*[dumbbell(:,2);dumbbell(1,2)],'c','Linewidth',2, 'LineSmoothing','on')

if (saveSimulations == true)
    saveas(gcf,fullfile(results_path,'OUDrift1.fig'))
    saveas(gcf,fullfile(results_path,'OUDrift1.png'))
end
pause(5)

view(3)
if (saveSimulations == true)
    saveas(gcf,fullfile(results_path,'OUDrift2.fig'))
    saveas(gcf,fullfile(results_path,'OUDrift2.png'))
end
pause(5)

% add plot of the origin dumbell to the second orientation



% --------------------------------- LA drift ------------------------------

% Notes: maybe make it equal to original so it gets preserved
% first run does not show the current state

s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s);


T = 30;
dt = 0.01;


%theta_LA = [0.5; 0.6];
%[bdryPts_t ctrlPts_t alpha_t] = Diffusion_drift_LA(x,x_small,800,100,theta_LA,10,10,dt,T);

theta_LA = [0.1; 0.05]; %[area_coef; length_coef]
[bdryPts_t ctrlPts_t alpha_t] = Diffusion_drift_LA(x,x_small,pi*100,2*pi*10,theta_LA,10,10,dt,T);

plot3D(bdryPts_t,ctrlPts_t,dt)
title('Shape Gradient Drift','FontSize',14,'fontweight','bold')

if (saveSimulations == true)
    saveas(gcf,fullfile(results_path,'ShapeDrift1.fig'))
    saveas(gcf,fullfile(results_path,'ShapeDrift1.png'))
end
pause(5)

view(3)
if (saveSimulations == true)
    saveas(gcf,fullfile(results_path,'ShapeDrift2.fig'))
    saveas(gcf,fullfile(results_path,'ShapeDrift2.png'))
end
pause(5)


% ---------------------------- regression-like OU -------------------------

s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s);


% Notes: not sure which parameters to put in (in notes dt = 0.1?)
% need to improve interpretation

T = 30;
dt = 0.01;

% ellipse_lr = 1.5*ellipse_lr;
% ellipse_ud = 1.5*ellipse_ud;

% [bdryPts_t ctrlPts_t alpha_t] = Diffusion_drift_OU_multiple(x,x_small, ellipse_lr,ellipse_ud,zeros(100,100),10,10,[0.05 0.05],dt,T);
% [bdryPts_t ctrlPts_t alpha_t] = Diffusion_drift_OU_multiple(x,x_small, ellipse_lr,ellipse_ud,zeros(100,100),10,10,[0.5 0.5],dt;,T);
[bdryPts_t ctrlPts_t alpha_t] = Diffusion_drift_OU_multiple(1.5*x,1.5*x_small, 1.5*ellipse_lr,1.5*ellipse_ud,zeros(100,100),10,10,[0.05 0.05],dt,T);

plot3D(bdryPts_t,ctrlPts_t,dt)
title('Regression-like Ornstein-Uhlenbeck Drift','FontSize',14, 'fontweight','bold')
plot3((T+1)*ones(size(ellipse_lr,1)+1,1),[1.5*ellipse_lr(:,1);1.5*ellipse_lr(1,1)],[1.5*ellipse_lr(:,2);1.5*ellipse_lr(1,2)],'c','LineWidth',2,'LineSmoothing','on')
plot3((T+1)*ones(size(ellipse_ud,1)+1,1),[1.5*ellipse_ud(:,1);1.5*ellipse_ud(1,1)],[1.5*ellipse_ud(:,2);1.5*ellipse_ud(1,2)],'c','LineWidth',2,'LineSmoothing','on')

if (saveSimulations == true)
    saveas(gcf,fullfile(results_path,'OUDrift_multiple1.fig'))
    saveas(gcf,fullfile(results_path,'OUDrift_multiple1.png'))
end
pause(5) 

view(3)
if (saveSimulations == true)
    saveas(gcf,fullfile(results_path,'OUDrift_multiple2.fig'))
    saveas(gcf,fullfile(results_path,'OUDrift_multiple2.png'))
end





% break



% -------------------------------------------------------------------
% Parameter Estimation
% -------------------------------------------------------------------

% reinitialize random stream

s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s);


% 1) Constant

theta_constant = 0.5*ones(size(x_small));
dt = 0.01;
T = 10;



[mse theta_hat] = MSE('constant',theta_constant,dt,T,100);

theta_hat_transformed = cat(4,theta_hat{:});

figure(1)

xaxis = linspace(dt,T,size(squeeze(theta_hat_transformed(1,1,:,:)),1))';
plot(xaxis,squeeze(theta_hat_transformed(1,1,:,:)),'LineSmoothing','on')

hold on
plot(xaxis,squeeze(theta_constant(1,1,:))*ones(size(xaxis)),'r-','LineWidth',2)
title('Constant Drift MLE (estimates for one coordinate)','FontSize',14, 'fontweight','bold')
xlabel('Time')
hold off

theta_constant
mean(theta_hat_transformed(:,:,end,:),4)

if (saveEstimations)
    saveas(gcf,fullfile(results_path,'constant_theta_hat.fig'))
    saveas(gcf,fullfile(results_path,'constant_theta_hat.png'))
end

% quantile plot
figure
cmap = colormap(cool(100));
set(gca, 'ColorOrder', cmap, 'NextPlot', 'replacechildren');
quantiles = [0.01,0.99, 0.05, 0.95,0.25,0.75];
bounds = abs(quantile(squeeze(theta_hat_transformed(1,1,:,:))',quantiles) - theta_constant(1,1));
[l,p] = boundedline(xaxis,repmat(theta_constant(1,1)*ones(size(xaxis))',length(quantiles)/2,1), reshape(bounds',[length(xaxis),2,length(quantiles)/2]),'cmap',cool(length(quantiles)/2));
outlinebounds(l,p);
title(sprintf('Constant Drift MLE Quantile (estimates for one coordinate)\n (0.01, 0.05, 0.25, 0.75, 0.95, 0.99)') ,'FontSize',14, 'fontweight','bold')
xlabel('Time')

if (saveEstimations)
    saveas(gcf,fullfile(results_path,'constant_theta_hat_qb.fig'))
    saveas(gcf,fullfile(results_path,'constant_theta_hat_qb.png'))
end



% 2) OU


s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s);

dt = 0.05;
T = 30;

[mse_OU theta_hat_OU] = MSE('OU',theta_OU,dt,T,100);

theta_hat_transformed = cat(3,theta_hat_OU{:});
figure(2)
cmap = colormap(cool(100));
set(gca, 'ColorOrder', cmap, 'NextPlot', 'replacechildren');
xaxis = linspace(dt,T,size(squeeze(theta_hat_transformed(1:end,:)),1))';
plot(xaxis,squeeze(theta_hat_transformed(:,:)),'LineSmoothing','on')

hold on
plot(xaxis,theta_OU*ones(size(xaxis)),'r-','LineWidth',2)
title('OU Drift MLE','FontSize',14, 'fontweight','bold')
xlabel('Time')
hold off

if (saveEstimations)
    saveas(gcf,fullfile(results_path,'OU_theta_hat.fig'))
    saveas(gcf,fullfile(results_path,'OU_theta_hat.png'))
end

disp(sprintf('The theta is %d.',theta_OU))
disp(sprintf('The estimated theta is %d.',mean(theta_hat_transformed(end,:,:),3)))


% quantile plot
figure
quantiles = [0.01,0.99, 0.05, 0.95,0.25,0.75];
bounds = abs(quantile(squeeze(theta_hat_transformed)',quantiles) - theta_OU);
[l,p] = boundedline(xaxis,repmat(theta_OU*ones(size(xaxis))',length(quantiles)/2,1), reshape(bounds',[length(xaxis),2,length(quantiles)/2]),'cmap',cool(length(quantiles)/2));
outlinebounds(l,p);
title(sprintf('OU Drift MLE Quantiles\n (0.01, 0.05, 0.25, 0.75, 0.95, 0.99)'),'FontSize',14, 'fontweight','bold')
xlabel('Time')

if (saveEstimations)
    saveas(gcf,fullfile(results_path,'OU_theta_hat_qb.fig'))
    saveas(gcf,fullfile(results_path,'OU_theta_hat_qb.png'))
end




% 3) LA
 
% Note: initial value is a huge negative or positive number:
% is the matrix signular?
% this is why I am not displaying the firs value in the plot (2:end)

s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s);

dt = 0.01;
T = 30;

[mse_LA theta_hat_LA] = MSE('LA',theta_LA,dt,T,100);

theta_hat_transformed = cat(4,theta_hat_LA{:});
figure(3)
hold off
cmap = colormap(cool(100));
set(gca, 'ColorOrder', cmap, 'NextPlot', 'replacechildren');
%subplot(2,1,1)
xaxis = linspace(dt,T,size(squeeze(theta_hat_transformed(1,2:end,:)),1))';
% plot(dt:dt:T-dt,squeeze(theta_hat_transformed(1,2:end,:)),'LineSmoothing','on')
plot(xaxis,squeeze(theta_hat_transformed(1,2:end,:)),'LineSmoothing','on')


hold on
plot(xaxis,theta_LA(2),'r.','LineWidth',2,'LineSmoothing','on')
title('Area Coefficient MLE','FontSize',14, 'fontweight','bold')
xlabel('Time')
hold off

if (saveEstimations)
    saveas(gcf,fullfile(results_path,'LA_theta_hat_area.fig'))
    saveas(gcf,fullfile(results_path,'LA_theta_hat_area.png'))
end
yl = ylim;


% quantile plot
figure
ylim(yl)

quantiles = [0.01,0.99, 0.05, 0.95,0.25,0.75];
bounds = abs(quantile(squeeze(theta_hat_transformed(1,2:end,:))',quantiles) - theta_LA(2));
[l,p] = boundedline(xaxis,repmat(theta_LA(2)*ones(size(xaxis))',length(quantiles)/2,1), reshape(bounds',[length(xaxis),2,length(quantiles)/2]),'cmap',cool(length(quantiles)/2));
outlinebounds(l,p);
title(sprintf('Area Coefficient MLE Quantiles\n  (0.01, 0.05, 0.25, 0.75, 0.95, 0.99)'),'FontSize',14, 'fontweight','bold')
xlabel('Time')

if (saveEstimations)
    saveas(gcf,fullfile(results_path,'LA_theta_hat_area_qb.fig'))
    saveas(gcf,fullfile(results_path,'LA_theta_hat_area_qb.png'))
end

figure(4)
%subplot(2,1,2)
hold off
ylim auto
cmap = colormap(cool(100));
set(gca, 'ColorOrder', cmap, 'NextPlot', 'replacechildren');

plot(xaxis,squeeze(theta_hat_transformed(2,2:end,:)),'LineSmoothing','on')

hold on
plot(xaxis,theta_LA(1)*ones(size(xaxis)),'r-','LineWidth',2)
title('Length Coefficient MLE','FontSize',14, 'fontweight','bold')
xlabel('Time')


if (saveEstimations)
    saveas(gcf,fullfile(results_path,'LA_theta_hat_length.fig'))
    saveas(gcf,fullfile(results_path,'LA_theta_hat_length.png'))
end

hold off

disp(sprintf('The theta is %d.',theta_LA(2)))
disp(sprintf('The estimated theta is %d.',mean(theta_hat_transformed(end,:,:),3)))


% quantile plot
figure
quantiles = [0.01,0.99, 0.05, 0.95,0.25,0.75];
bounds = abs(quantile(squeeze(theta_hat_transformed(2,2:end,:))',quantiles) - theta_LA(1));
[l,p] = boundedline(xaxis,repmat(theta_LA(1)*ones(size(xaxis))',length(quantiles)/2,1), reshape(bounds',[length(xaxis),2,length(quantiles)/2]),'cmap',cool(length(quantiles)/2));
outlinebounds(l,p);
title(sprintf('Length Coefficient MLE Quantiles \n (0.01, 0.05, 0.25, 0.75, 0.95, 0.99)'),'FontSize',14, 'fontweight','bold')
xlabel('Time')

if (saveEstimations)
    saveas(gcf,fullfile(results_path,'LA_theta_hat_length_qb.fig'))
    saveas(gcf,fullfile(results_path,'LA_theta_hat_length_qb.png'))
end





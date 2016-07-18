% This script stores the experiments for simulating diffusions and estimating their parameters

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

[bdryPts_t, ctrlPts_t] = Diffusion(x,x_small,zeros(size(x_small)),10,0.05,30);
plot3D(bdryPts_t,ctrlPts_t,0.05)

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

theta_constant = 0.1*randn(size(x_small));
theta_constant = 0.5*ones(size(x_small));
[bdryPts_t, ctrlPts_t] = Diffusion(x,x_small,theta_constant,10,0.001,30);

plot3D(bdryPts_t,ctrlPts_t,0.001)
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
dt = 0.05; % not sure of 0.05 or 0.01 is better
T = 30;
theta_OU = 0.5;
[bdryPts_t ctrlPts_t alpha_t] = Diffusion_drift_OU(1.5*ellipse_lr,1.5*ellipse_lr(1:5:end-4,:),1.5*dumbbell,zeros(100,100),10,10,theta_OU,dt,T);

plot3D(bdryPts_t,ctrlPts_t,dt)
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


%theta_LA = [0.5; 0.6];
%[bdryPts_t ctrlPts_t alpha_t] = Diffusion_drift_LA(x,x_small,800,100,theta_LA,10,10,dt,T);

theta_LA = [0.1; 0.05]; %[area_coef; length_coef]
[bdryPts_t ctrlPts_t alpha_t] = Diffusion_drift_LA(x,x_small,pi*100,2*pi*10,theta_LA,10,10,0.01,30);

plot3D(bdryPts_t,ctrlPts_t,0.01)
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

% Notes: not sure which parameters to put in (in notes dt = 0.1?)
% need to improve interpretation

[bdryPts_t ctrlPts_t alpha_t] = Diffusion_drift_OU_multiple(x,x_small, ellipse_lr,ellipse_ud,zeros(100,100),10,10,[0.05 0.05],0.01,T);

plot3D(bdryPts_t,ctrlPts_t,0.01)
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
pause(5)




break



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
plot(dt:dt:T,squeeze(theta_hat_transformed(1,1,:,:)))

hold on
plot(dt:dt:T,squeeze(theta_constant(1,1,:)),'r.','LineWidth',3)
title('Constant Drift MSE (estimates for one coordinate)')
hold off

theta_constant
mean(theta_hat_transformed(:,:,end,:),4)

if (saveEstimations)
    saveas(gcf,fullfile(results_path,'constant_theta_hat.fig'))
    saveas(gcf,fullfile(results_path,'constant_theta_hat.png'))
end



% 2) OU

dt = 0.05;
T = 30;

[mse_OU theta_hat_OU] = MSE('OU',theta_OU,dt,T,100);

theta_hat_transformed = cat(3,theta_hat_OU{:});
figure(2)
plot(dt:dt:T,squeeze(theta_hat_transformed(:,:)))

hold on
plot(dt:dt:T,theta_OU,'r.','LineWidth',3)
title('OU Drift')
hold off

if (saveEstimations)
    saveas(gcf,fullfile(results_path,'OU_theta_hat.fig'))
    saveas(gcf,fullfile(results_path,'OU_theta_hat.png'))
end

disp(sprintf('The theta is %d.',theta_OU))
disp(sprintf('The estimated theta is %d.',mean(theta_hat_transformed(end,:,:),3)))



% 3) LA
 
% Note: initial value is a huge negative or positive number:
% is the matrix signular?
% this is why I am not displaying the firs value in the plot (2:end)

dt = 0.01;
T = 30;

[mse_LA theta_hat_LA] = MSE('LA',theta_LA,dt,T,100);

theta_hat_transformed = cat(4,theta_hat_LA{:});
figure(3)
%subplot(2,1,1)
plot(dt:dt:T,squeeze(theta_hat_transformed(2,2:end,:)))

hold on
plot(dt:dt:T,theta_LA(1),'r.','LineWidth',3)
title('Area Coefficient')
hold off

if (saveEstimations)
    saveas(gcf,fullfile(results_path,'LA_theta_hat_area.fig'))
    saveas(gcf,fullfile(results_path,'LA_theta_hat_area.png'))
end

figure(4)
%subplot(2,1,2)

plot(dt:dt:T,squeeze(theta_hat_transformed(2,2:end,:)))

hold on
plot(dt:dt:T,theta_LA(2),'r.','LineWidth',3)
title('Length Coefficient')
hold off

if (saveEstimations)
    saveas(gcf,fullfile(results_path,'LA_theta_hat_length.fig'))
    saveas(gcf,fullfile(results_path,'LA_theta_hat_length.png'))
end

disp(sprintf('The theta is %d.',theta_OU))
disp(sprintf('The estimated theta is %d.',mean(theta_hat_transformed(end,:,:),3)))





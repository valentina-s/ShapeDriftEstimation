% This script stores the experiments for simulating diffusions and estimating their parameters

s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s);

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
pause(5)
view(3)
pause(5)


% ----------------------------- constant drift ----------------------------
% Notes: make it more drifty

theta_constant = 0.1*randn(size(x_small));
[bdryPts_t, ctrlPts_t] = Diffusion(x,x_small,theta_constant,10,0.05,30);
plot3D(bdryPts_t,ctrlPts_t,0.05)
pause(5)
view(3)
pause(5)

% -------------------------------- OU drift -------------------------------
dt = 0.05;
T = 30;
theta_OU = 0.5;
[bdryPts_t ctrlPts_t alpha_t] = Diffusion_drift_OU(1.5*ellipse_lr,1.5*ellipse_lr(1:5:end-4,:),1.5*dumbbell,zeros(100,100),10,10,theta_OU,dt,T);
plot3D(bdryPts_t,ctrlPts_t,0.05)
pause(5)
view(3)
pause(5)



% --------------------------------- LA drift ------------------------------

% Notes: maybe make it equal to original so it gets preserved
% first run does not show the current state


%theta_LA = [0.5; 0.6];
%[bdryPts_t ctrlPts_t alpha_t] = Diffusion_drift_LA(x,x_small,800,100,theta_LA,10,10,dt,T);

theta_LA = [0.1; 0.05]; %[area_coef; length_coef]
[bdryPts_t ctrlPts_t alpha_t] = Diffusion_drift_LA(x,x_small,pi*100,2*pi*10,theta_LA,10,10,0.01,30);
plot3D(bdryPts_t,ctrlPts_t,0.01)
pause(5)
view(3)
pause(5)


% ---------------------------- regression-like OU -------------------------

% Notes: not sure which parameters to put in (in notes dt = 0.1?)
% need to improve interpretation

[bdryPts_t ctrlPts_t alpha_t] = Diffusion_drift_OU_multiple(x,x_small, ellipse_lr,ellipse_ud,zeros(100,100),10,10,[0.05 0.05],0.01,T);
plot3D(bdryPts_t,ctrlPts_t,0.01)
pause(5) 
view(3)
pause(5)








% -------------------------------------------------------------------
% Parameter Estimation
% -------------------------------------------------------------------

% 1) Constant

[mse theta_hat] = MSE('constant',theta_constant,0.05,30,100);

theta_hat_transformed = cat(4,theta_hat{:});

figure(1)
plot(squeeze(theta_hat_transformed(1,1,:,:)))

hold on
plot(1:600,squeeze(theta_constant(1,1,:)),'r.','LineWidth',3)
title('Constant Drift MSE (estimates for one coordinate)')
hold off



% 2) OU

[mse_OU theta_hat_OU] = MSE('OU',theta_OU,0.05,30,100);

theta_hat_transformed = cat(3,theta_hat_OU{:});
figure(2)
plot(squeeze(theta_hat_transformed(:,:)))

hold on
plot(1:600,theta_OU,'r.','LineWidth',3)
title('OU Drift')
hold off


% 3) LA
 
% Note: initial value is a huge negative or positive number:
% is the matrix signular?
% this is why I am not displaying the firs value in the plot (2:end)

[mse_LA theta_hat_LA] = MSE('LA',theta_LA,0.01,30,100);

theta_hat_transformed = cat(4,theta_hat_LA{:});
figure(3)
subplot(2,1,1)
plot(squeeze(theta_hat_transformed(2,2:end,:)))

hold on
plot(1:size(theta_LA,1),theta_LA(1),'r.','LineWidth',3)
title('Area Coefficient')
hold off

subplot(2,1,2)

plot(squeeze(theta_hat_transformed(2,2:end,:)))

hold on
plot(1:size(theta_LA,1),theta_LA(2),'r.','LineWidth',3)
title('Length Coefficient')
hold off




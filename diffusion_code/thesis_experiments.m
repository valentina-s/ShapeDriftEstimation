% This script stores the experiments for simulating diffusions and estimating their parameters

s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s);

% Diffusion Simulation

% driftless

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

% no drift
[bdryPts_t, ctrlPts_t] = Diffusion(x,x_small,zeros(size(x_small)),10,0.05,30);
plot3D(bdryPts_t,ctrlPts_t,0.05)
pause(5)
view(3)
pause(5)

% Notes: make it more drifty
% constant drift
[bdryPts_t, ctrlPts_t] = Diffusion(x,x_small,0.1*randn(size(x_small)),10,0.05,30);
plot3D(bdryPts_t,ctrlPts_t,0.05)
pause(5)
view(3)
pause(5)

% OU drift
dt = 0.05;
T = 30;
[bdryPts_t ctrlPts_t alpha_t] = Diffusion_drift_OU(1.5*ellipse_lr,1.5*ellipse_lr(1:5:end-4,:),1.5*dumbbell,zeros(100,100),10,10,0.5,dt,T);
plot3D(bdryPts_t,ctrlPts_t,0.05)
pause(5)
view(3)
pause(5)


% Notes: maybe make it equal to original so it gets preserved
% first run does not show the current state
% LA drift
[bdryPts_t ctrlPts_t alpha_t] = Diffusion_drift_LA(x,x_small,800,100,0.01,0.1,10,10,dt,T);
plot3D(bdryPts_t,ctrlPts_t,0.05)
pause(5)
view(3)
pause(5)

% Notes: not sure which parameters to put it (dt = 0.1?), not working as expected
% regression-like OU
[bdryPts_t ctrlPts_t alpha_t] = Diffusion_drift_OU_multiple(x,x_small, ellipse_lr,ellipse_ud,zeros(100,100),10,10,[0.5 0.5],dt,T);
plot3D(bdryPts_t,ctrlPts_t,0.05)
pause(5)
view(3)
pause(5)









% Parameter Estimation

% 1) Constant

% [mse theta_hat] = MSE('constant',theta,0.05,10,4);

% 2) OU

% [mse_OU theta_hat_OU] = MSE('OU',2,0.05,30,4);


% 3) LA

% [mse_LA theta_hat_LA] = MSE('LA',[0.5,0.6],0.05,30,4);


function [mse theta_hat] = MSE(drift_type,theta,dt,T,sampleSize)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% [mse theta_hat] = MSE(drift_type,theta,dt,T,sampleSize)
% 
% calculate MSE of the likelihood ratio estimator
% 
% INPUTS:
%   drift_type - a string indicating the type of the drift
%               'constant', 'OU', 'LA'
%   theta      - the coefficients with which to simulate the diffusions
%                1D for 'constant'
%                1D for 'OU'
%                2D for 'LA'
%   dt         - time step
%   T          - final time
%   sampleSize - number of paths
%
%
% OUTPUTS:
%   mse        - mse value
%   theta_hat  - the estimated values of theta - cell of size sampleSize
%                each entry containing the theta estimate at each time
%                sizeOfTheta x # time steps
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initializing variables
total = 0;
s = 0;
err = zeros(sampleSize,1);

% the dimension of theta depends on the type of drift
% theta_hat = zeros(size(theta,1),size(theta,2),sampleSize);
clear theta_hat

a = [0:0.1:2*pi];
x = 10*[cos(a)' sin(a)'] + 30;
x_small = x(1:7:end,:);
circle = [5*sin(a)' 5*cos(a)'] + 30;
ellipse_lr = [10*sin(a)' 5*cos(a)'] + 30;
ellipse_ud = [5*sin(a)' 10*cos(a)'] + 30;

a = 0.9;
y = [-1:0.05:1];
yy_plus = (y.^4.*(a - y.^2)/(a^2));
yy_minus = -y.^4.*(a - y.^2)/(a^2);
dumbbell = [[y';y(end:-1:1)'] [yy_plus'+0.15 ; yy_minus'-0.15]]*15+30;


%  switch drift_type
%      case 'constant'
%          theta_hat = zeros(size(x_small,1),2,sampleSize);
%      case 'OU'
%          theta_hat = zeros(sampleSize,1);
%      case 'LA'
%          theta_hat = zeros(2,sampleSize);
%  end

theta_hat = cell(sampleSize,1);



% generate sample paths
parfor k=1:sampleSize
 

    disp(sprintf('Iteration number %i',k))
    
    switch drift_type    
  
    case 'constant'
        [bdryPts_t ctrlPts_t alpha_t] = Diffusion(x,x_small,theta,10,dt,T);
        theta_hat{k} = drift_estimate(ctrlPts_t,alpha_t,dt,10);
        err(k,1) = sum(sum((theta_hat{k}(:,:,end) - theta).^2));
        

    case 'OU'
        [bdryPts_t ctrlPts_t alpha_t] = Diffusion_drift_OU(1.5*ellipse_lr,1.5*ellipse_lr(1:5:end-4,:),1.5*dumbbell,zeros(100,100),10,10,theta,dt,T);
        theta_hat{k} = drift_OU_estimate(bdryPts_t,ctrlPts_t,1.5*dumbbell,alpha_t,dt,10);
        err(k,1) = sum(sum((theta_hat{k}(end) - theta).^2));
  
    case 'LA'
        [bdryPts_t ctrlPts_t alpha_t] = Diffusion_drift_LA(x,x_small,800,100,theta,10,10,dt,T);
        theta_hat{k} = drift_LA_estimate(bdryPts_t,ctrlPts_t,800,100,alpha_t,dt,10);
        disp(size(theta_hat{k}))
        err(k,1) = sum(sum((theta_hat{k}(:,end) - theta).^2));
    end    
%     
%  end
%  final_theta_hat = sum(theta_hat,3)/sampleSize;

end
mse = sum(err)/sampleSize;




    

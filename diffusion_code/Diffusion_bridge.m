function [x_t alpha_t] = Diffusion_bridge(x0,x1,alpha_drift,sigma,dt,T)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% [x_t alpha_t] = Diffusion_bridge(x0,x1,alpha_drift,sigma,dt,T)
%
% x0 - the initial configuration of the landmarks
% drift - contains the coefficients of the
% sigma - the deformation kernel width
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% this parameter determines 
plotting = false;

N = T/dt;
x = x0;
nofCtrlPts = size(x0,1);
x_t = zeros(nofCtrlPts,2,N);
x_t(:,:,1) = x0;
alpha_t = zeros(nofCtrlPts,2,N);


if plotting
% figure(1)
  hold off
  h1 = plot(x0(:,1),x0(:,2),'co--')
  plot(x1(:,1),x1(:,2),'bo--')
  hold on
end

for i = 1:N-1
    % generate a random tangent vector at x
    K = K_matrix(x_t(:,:,i),x_t(:,:,i),sigma*sigma);
    % alpha_t(:,:,i) = inv(chol(K))*alpha_drift + (1/sqrt(dt))*inv(chol(K))*randn(size(alpha_drift));
    % move along the exponential map
    % x_t(:,:,i+1) = exponentialR(x_t(:,:,i),x_t(:,:,i),alpha_t(:,:,i),sigma,dt,10);
    %K = eye(nofCtrlPts);
    bridge = (x1 - x_t(:,:,i))/(T-i*dt);
    
    x_t(:,:,i+1) = x_t(:,:,i) + dt*chol(K)*alpha_drift + dt*bridge + sqrt(dt) * sqrtm(K)*randn(size(alpha_drift));
    
    if plotting

        set(h1,'Color','b')

        h1 = plot(x_t(:,1,i+1),x_t(:,2,i+1),'ro--', 'Linewidth',2)
        pause(0.01)

    end

end
 plot(x_t(:,1,i+1),x_t(:,2,i+1),'go--')



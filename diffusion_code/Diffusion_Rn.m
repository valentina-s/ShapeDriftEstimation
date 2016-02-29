function [x_t alpha_t] = Diffusion_Rn(x0,alpha_drift,sigma,dt,T)
% x0 - the initial configuration of the landmarks
% drift - contains the coefficients of the
% sigma - the deformation kernel width

plotting = true;

N = T/dt;
x = x0;
nofCtrlPts = size(x0,1);
x_t = zeros(nofCtrlPts,2,N+1);
x_t(:,:,1) = x0;
alpha_t = zeros(nofCtrlPts,2,N);



figure(1)
hold off
plot(x(:,1),x(:,2),'ro--')
hold on

for i = 1:N
    % generate a random tangent vector at x
    K = K_matrix(x_t(:,:,i),x_t(:,:,i),sigma*sigma);
    %K = eye(size(K));
    % alpha_t(:,:,i) = inv(chol(K))*alpha_drift + (1/sqrt(dt))*inv(chol(K))*randn(size(alpha_drift));
    % move along the exponential map
    % x_t(:,:,i+1) = exponentialR(x_t(:,:,i),x_t(:,:,i),alpha_t(:,:,i),sigma,dt,10);
    %K = eye(nofCtrlPts);
    x_t(:,:,i+1) = x_t(:,:,i) + dt*sqrtm(K)*alpha_drift + sqrt(dt) * (sqrt(K))*randn(size(alpha_drift));
    
    if plotting 
      plot(x_t(:,1,i+1),x_t(:,2,i+1),'ro')
      hold on
      pause(0.01)
    end
end

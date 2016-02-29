function [theta_hat total] = drift_estimate(x_t,alpha_t,dt,sigma)

% x_t is the sequence of observations 

% rough version - dX = X_{t+1} - X_{t}
if isempty(alpha_t)
    dx_t = circshift(x_t,[0 0 -1]) - x_t;
else
    dx_t = alpha_t;
end
N = size(dx_t,3);
total = 0;

% calculate the average
theta_hat = zeros(size(dx_t,1),2,N-1);
for i = 1:N-1
    if  isempty(alpha_t)
        K = chol(K_matrix(x_t(:,:,i),x_t(:,:,i),sigma^2));
        %K = eye(size(dx_t,1));
        total = total + K'\(x_t(:,:,i+1) - x_t(:,:,i));
    
        theta_hat(:,:,i) = total/(dt*i);
    else
    
    
        K = chol(K_matrix(x_t(:,:,i),x_t(:,:,i),sigma^2));
        %K = eye(size(K));
        total = total + K*dx_t(:,:,i);
        theta_hat(:,:,i) = total/(dt*i);
    end
end



% equivalently
% theta_hat = (x_t(:,:,end) - x_t(:,:,1))/(dt*(N-1));





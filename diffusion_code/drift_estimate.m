function [theta_hat total] = drift_estimate(x_t,alpha_t,dt,sigma)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% [theta_hat total] = drift_estimate(x_t,alpha_t,dt,sigma)
% 
% INPUTS:
%   x_t is the sequence of observations
%   alpha_t
%   dt
%   sigma
%
% OUTPUTS:
%   theta_hat
%   total
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If alpha_t is not provided it is estimated
% rough version - dX = X_{t+1} - X_{t}
if isempty(alpha_t)
    dx_t = circshift(x_t,[0 0 -1]) - x_t;
    N = size(dx_t,3);
    d = size(dx_t,1);
else 
    N = size(alpha_t,3);
    d = size(alpha_t,1);
end

total = 0;

% calculate the average
theta_hat = zeros(d,2,N-1);
for i = 1:N-1
    if  isempty(alpha_t)
        K = K_matrix(x_t(:,:,i),x_t(:,:,i),sigma^2);
        % K = eye(d);
        total = total + K'\(x_t(:,:,i+1) - x_t(:,:,i));
        cond(K)
        
        %total = total + inv(K)*(x_t(:,:,i+1) - x_t(:,:,i));
    
        theta_hat(:,:,i) = total/(dt*i);
    else
    
    
        % K = chol(K_matrix(x_t(:,:,i),x_t(:,:,i),sigma^2));
        % K = eye(d);
        % total = total + K*alpha_t(:,:,i);
        total = total + alpha_t(:,:,i);
        theta_hat(:,:,i) = total/(dt*i);
    end
end



% equivalently
% theta_hat = (x_t(:,:,end) - x_t(:,:,1))/(dt*(N-1));





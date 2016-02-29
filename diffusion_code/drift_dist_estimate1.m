function theta_hat = drift_dist_estimate1(x_t,xx_t,x_mean,alpha_t,dt,sigma)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function theta_hat = drift_OU_estimate(x_t,xx_t,x_mean,alpha_t,dt,sigma)
%
%
% x_t is the sequence of observations 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% rough version - dX = X_{t+1} - X_{t}

dx_t = circshift(x_t,[0 0 -1]) - x_t;

N = size(dx_t,3);


% calculate the average

% store the matrix in M and the vector in b

M = 0;
b = 0;
theta_hat =0;
for i = 1:N-1

    % calculate drift term with respect to the distance function
    g =  x_mean - x_t(:,:,i);
    
    M = M + g(:)'*g(:)*dt;

    temp = dx_t(:,:,i);

    b = b + g(:)'*temp(:);

end

theta_hat =  b/M;

end


function theta_hat = drift_dist_estimate(x_t,xx_t,x_mean,alpha_t,dt,sigma)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function theta_hat = drift_OU_estimate(x_t,xx_t,x_mean,alpha_t,dt,sigma)
%
%
% x_t is the sequence of observations 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% rough version - dX = X_{t+1} - X_{t}
if isempty(alpha_t)
    dx_t = circshift(x_t,[0 0 -1]) - x_t;
else
    dx_t = alpha_t;
end
N = size(dx_t,3);


% calculate the average

% store the matrix in M and the vector in b

M = 0;
b = 0;
theta_hat =0;
for i = 1:N-1

    % calculating the kernel matrices
    K_bdry = K_matrix(x_t(:,:,i),xx_t(:,:,i),sigma^2);
    K_ctrl = K_matrix(xx_t(:,:,i),xx_t(:,:,i),sigma^2);


    % calculate drift term with respect to the distance function
    g =  x_mean - x_t(:,:,i);
    g1  = K_bdry'*g; 


    galpha(:,1) = K_ctrl\g1(:,1);
    galpha(:,2) = K_ctrl\g1(:,2);

    drift = galpha(:);
    
    M = M + drift'*blkdiag(inv(K_ctrl),inv(K_ctrl))*drift*dt;

    if isempty(alpha_t)
       alpha = K_bdry'\dx_t(:,:,i);
       temp = K_ctrl*alpha;
    else
       temp = K_ctrl*alpha_t(:,:,i);
    end


    b = b + drift'*blkdiag(inv(K_ctrl),inv(K_ctrl))*temp(:);

    %if  isempty(alpha_t)
    %    K = chol(K_matrix(xx_t(:,:,i),xx_t(:,:,i),sigma^2));
    %    
    %     + K'\(x_t(:,:,i+1) - x_t(:,:,i));
    %else
    %    K = chol(K_matrix(xx_t(:,:,i),xx_t(:,:,i),sigma^2));
    % end
end

theta_hat = - b/M;

end



function [x_t] = OU_bridge(x0,x1,mu,theta,dt,T,sampleSize)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function [x_t] = OU_bridge(x0,x1,mu,theta,dt,T,sampleSize)
%
% OU_bridge generates a sample from the bridge of the OU process.
% 
% Note: it does not return the final destination point.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% number of observations
N = round(T/dt);
x_t = zeros(size(mu,1),N,sampleSize);
for k = 1:sampleSize
  x_t(:,1,k) = x0;
  for t = 1:N-1
    bridge = theta*(x1 - mu - (x_t(:,t,k)-mu)*exp(-theta*(T - t*dt)))*exp(-theta*(T-t*dt))/(1-exp(-2*theta*(T-t*dt)));
    x_t(:,t+1,k) = x_t(:,t,k) + theta*(mu - x_t(:,t,k))*dt + bridge*dt + sqrt(dt)*randn(size(mu));
  end
end

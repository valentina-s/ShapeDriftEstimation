function [x_t t] = OU_R(x0,mu,theta,dt,T)
N = round(T/dt); % modify if not an integer
t = [0:dt:T];
x_t = zeros(size(mu,1),N);
x_t(:,1) = x0;
for k =1:N
  x_t(:,k+1) = x_t(:,k) + theta*dt*(mu - x_t(:,k)) + sqrt(dt)*randn(size(mu,1),1);
end

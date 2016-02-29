N = 10;
MSE = zeros(N,1);
for n = 1:N
     [x_t t_axis] = OU_R(x0,ones(100,1),1,0.01,10);
     theta_hat = EM_OU_parallel(x_t(:,1:100:end),t_axis(1:100:end),ones(100,1),0.01,1,100,20,'exact');
     MSE(n,1) = theta_hat(end); 
end

function theta_hat = brownian_test(sampleSize,N)


drift = 1;
T = 10;
dt_small = 0.01;
nofIt = 50;

theta = 5;
x_t = diffusion(drift,theta,dt_small,T);

plot(x_t,'r')
hold on

x_sparse = x_t(1:N:end);
nofObs = size(x_sparse);

theta_hat = zeros(nofIt,1);
theta_hat(1) = 0.5;
for it = 1:nofIt

logL = zeros(sampleSize,1);
nom = zeros(sampleSize,1);
denom = zeros(sampleSize,1);

x_dense = zeros(size(x_t,1),sampleSize);

for k=1:sampleSize
    for t=1:size(x_sparse,1)-1   

        x_dense((t-1)*N+1:t*N,k) = brownianBridge(x_sparse(t),x_sparse(t+1),dt_small,N*dt_small);


    end

    x_dense(end,k) = x_sparse(end);



    logL(k) = logLikelihood(x_dense(:,k), drift, theta_hat(it), dt_small);
% logL(k)
% plot(x_dense(:,k))
% hold on
% pause(1)
    [dummy nom(k) denom(k)] = theta_estimate(x_dense(:,k),drift,dt_small);
    
   
end
nom
denom
logL



 L = exp(logL-max(logL));

 theta_hat(it+1) = (nom'*L)/(denom'*L);

end

% simulating brownian motion with drift
function x_t = diffusion(drift,theta,dt,T)

N = T/dt;
x_t = zeros(N+1,1);
for i=1:N
    x_t(i+1) = x_t(i) + dt*drift*theta + sqrt(dt)*randn(1);
end


% sampling brownian bridge
function x_t = brownianBridge(x0,x1, dt, T)
    N = T/dt;
    x_t = zeros(N,1);
    x_t(1) = x0;
    for i=1:N-1
        x_t(i+1) = x_t(i) + dt*(x1 - x_t(i))/(T - dt*(i)) + sqrt(dt)*randn(1);
    end


% estimating theta
function [theta, nom, denom] = theta_estimate(x_t,drift,dt)
    N = size(x_t,1);
    dx_t = circshift(x_t,[-1]) - x_t;
  
    nom = 0;
    denom = 0;

    for k = 1:N-1
        v = dx_t(k);
       

        % calculate the drift
        nom = nom + drift*v;
        denom = denom + drift*drift*dt; 

    end
    theta = nom/denom;

% calculating the likelihood ratio
function logL = logLikelihood(x_t, drift, theta,dt)
    logL = 0;
    



    N = size(x_t,1);
    dx_t = circshift(x_t,[-1]) - x_t;
    T = N*dt;


for k = 2:N
  v = dx_t(k-1);
  x = x_t(k);
  x_final = x_t(end);
 
  type = 'drift';
  
  switch type
  case 'drift'
   % proposal with drift 
   logL = logL  + (x - x_final)*drift/(N*dt - (k-1)*dt)*dt;
  case 'driftless'
    % proposal without drift
    logL = logL + theta*drift*v  - theta*theta*drift*drift*dt/2;
  end
end
    
    
    

    
    






    

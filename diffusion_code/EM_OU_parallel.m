function theta_hat = EM_OU_parallel(bdryPts_t,t_axis,mu,dt_small,dt_large,sampleSize,nofIt,type)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% theta_hat = EM_OU(bdryPts_t,mu,dt_small,dt_large,sampleSize,nofIt)
%
% This program implements an EM algorithm to estimate a parameter theta
% from a discrete sequence of observations. 
%
% bdryPts_t -	 the discrete sequence of observations
% dt_small -     the time step between every two dense observations
% dt_large -     the time step between every two sparse observations
% T -		 the distance between every two discrete observations
% sampleSize -	 the number of samples in the approximation of the expectation
% nofIt -	 number of EM iterations
% type - type of the proposal - 'exact','driftless','drift'
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nofBdryPts = size(bdryPts_t,1);
% the number of dense observations to generate between any two discrete observations
N = round(dt_large/dt_small); 
nofObs = size(bdryPts_t,2);
theta_hat = zeros(nofIt+1,1);
x_dense = zeros(nofBdryPts,nofObs*N);
alpha_dense = zeros(nofBdryPts,nofObs*N);

% obtaining initial estimate for theta based on the sparse observations
theta_hat(1,1) = theta_estimate_OU(bdryPts_t,mu,dt_large);

sprintf('The initial estimate of theta is %d',theta_hat(1,1))


temp = repmat(bdryPts_t,[1 1 sampleSize]);



for it = 1:nofIt 
  disp(sprintf('Iteration %i',it))
  total = 0;
  s = zeros(nofBdryPts,2,sampleSize);
  
  logL = zeros(sampleSize,	1);  
  nom = zeros(sampleSize,1);
  denom = zeros(sampleSize,1);
 
tic 
  parfor k = 1:sampleSize  
    % disp(sprintf('Sample %i',k))
    % simulating the proposal paths between the discrete observations
    % size(sampleBridges(temp(:,:,k),mu,dt_small,dt_large,type))
 

	[logL(k) nom(k) denom(k)] = sampleBridges(temp(:,:,k),t_axis,mu,dt_small,dt_large,theta_hat(it,1),type);
	
	
  end
toc
  

  if strcmp(type,'drift') || strcmp(type,'driftless')
        L = exp(logL - max(logL));
        theta_hat(it+1) = nom'*L/(denom'*L);
  elseif strcmp(type,'exact')
        theta_hat(it+1) = sum(nom)/sum(denom);
  else
        error('Wrong proposal type.')
  end
  disp(sprintf('The estimate for theta is %d',theta_hat(it+1)))


  
end
end


function logL = logLikelihood(x_t,theta,mu,dt,type)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  logL = logLikelihood(x_t,theta,mu,dt,type)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


N = size(x_t,2);
dx_t = circshift(x_t,[0 -1]) - x_t;
T = N*dt;

% calculating the likelihood along the path
logL = 0;
for k = 2:N
  v = dx_t(:,k-1);
  x = x_t(:,k);
  % is this really the final value?
  x_final = x_t(:,end);
  drift = theta*(x_t(:,k) - mu);
  
  % calculate the drift
  term1 = drift(:)'*v(:);
  term2 = drift(:)'*drift(:);
  term3 = (x(:) - x_final(:))'*(x(:)-x_final(:))/(2*(T - dt*(k-1)));
  term4 = (x(:) - x_final(:))'*drift(:)/(T - dt*(k-1));
  % dont know if the third term is -/+ 1/2?
  
  switch type
  case 'drift'
   % proposal with drift 
   logL = logL  - term3 - term4*dt;
  case 'driftless'
    % proposal without drift
    logL = logL + term1 - term2*dt/2 - term3;
  end
   %logL
end

end


function logL = logLikelihood1(x_t,t_axis,theta,mu,dt,type)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  logL = logLikelihood(x_t,theta,mu,dt,type)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


N = size(x_t,2);
dx_t = circshift(x_t,[0 -1]) - x_t;
T = N*dt;

% calculating the likelihood along the path
logL = 0;
for k = 2:N
  v = dx_t(:,k-1);
  x = x_t(:,k);
  % is this really the final value?
  x_final = x_t(:,end);
  drift = theta*(x_t(:,k) - mu);
  
  % calculate the drift
  term1 = drift(:)'*v(:);
  term2 = drift(:)'*drift(:);
  term3 = (x(:) - x_final(:))'*(x(:)-x_final(:))/(2*(t_axis(end) - t_axis(k-1)));
  term4 = (x(:) - x_final(:))'*drift(:)/(t_axis(end) - t_axis(k-1));
  % dont know if the third term is -/+ 1/2?
  
  switch type
  case 'drift'
   % proposal with drift 
   logL = logL  + term4*dt;
  case 'driftless'
    % proposal without drift
    logL = logL + term1 - term2*dt/2;
  end
   %logL
end

end


function [theta, nom, denom] = theta_estimate(x_t,mu,dt)
N = size(x_t,2);
dx_t = circshift(x_t,[0 -1]) - x_t;
%  disp(dx_t(:,1))
  
nom = 0;
denom = 0;

for k = 2:N
  v = dx_t(:,k-1);
  drift = (mu - x_t(:,k-1));
  
  % calculate the drift
  nom = nom + drift(:)'*v(:);
  denom = denom + drift(:)'*drift(:)*dt; 

end
theta = nom/denom;
end


function [logL nom denom] = sampleBridges(bdryPts_t,t_axis,mu,dt_small,dt_large,theta_hat,type)
nofObs = size(bdryPts_t,2);
nofBdryPts = size(bdryPts_t,1);
N = dt_large/dt_small;
x_dense = zeros(nofBdryPts,size(0:dt_small:t_axis(end),2));



logL = 0;
x_dense(:,1) = bdryPts_t(:,1);
for t = 1:nofObs-1
                
        % sampling from proposal distribution
        switch type

        case 'driftless'
	        for j=1:N-1   
	            x_dense(:,(t-1)*N+j+1) = x_dense(:,(t-1)*N+j) + ...
                dt_small* (bdryPts_t(:,t+1) - x_dense(:,(t-1)*N+j))/(N*dt_small - j*dt_small) + ...
                sqrt(dt_small)*randn(size(mu));
	        end

	    case 'drift'
	        for j=1:N-1
	            drift = theta_hat*(mu - x_dense(:,(t-1)*N+j));
	            x_dense(:,(t-1)*N+j+1) = x_dense(:,(t-1)*N+j) + dt_small*drift + ...
                dt_small* (bdryPts_t(:,t+1) - x_dense(:,(t-1)*N+j))/(N*dt_small - j*dt_small) + ...
                sqrt(dt_small)*randn(size(mu));
	        end


        case 'exact'
            x_dense(:,(t-1)*N+1:t*N) = OU_bridge(bdryPts_t(:,t),bdryPts_t(:,t+1),mu,theta_hat,dt_small,dt_large,1);

        x_dense(1,t*N)
        bdryPts_t(1,t+1)
        end
        x_dense(:,t*N+1) = bdryPts_t(:,t+1);
	     plot(0:dt_small:t_axis(end),x_dense(1,:))
         hold on
         plot(t_axis,bdryPts_t(1,:),'ro','Linewidth',2)
                
         pause(0.1)
         hold off

	    % logL = logL + logLikelihood(x_dense(:,(t-1)*N+1:t*N+1),theta_hat,mu,dt_small,type);
       

end

logL = logLikelihood1(x_dense,0:dt_small:t_axis(end),theta_hat,mu,dt_small,type);

[dummy, nom, denom] = theta_estimate(x_dense,mu,dt_small);
end


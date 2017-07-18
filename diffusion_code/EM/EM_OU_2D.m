function theta_hat = EM_OU_2D(bdryPts_t,mu,sigma,dt_small,dt_large,sampleSize,nofIt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% theta_hat = EM_OU(bdryPts_t,mu,dt_small,dt_large,sampleSize,nofIt)
%
% This program implements an EM algorithm to estimate a parameter theta
% from a discrete sequence of observations. 
%
% bdryPts_t -	 the discrete sequence of observations
% dt -		 the distance between every two dense observations
% T -		 the distance between every two discrete observations
% sampleSize -	 the number of samples in the approximation of the expectation
% nofIt -	 number of EM iterations
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
type = 'drift';
nofBdryPts = size(bdryPts_t,1);
% the number of dense observations to generate between any two discrete observations
N = round(dt_large/dt_small); 
nofObs = size(bdryPts_t,3);
theta_hat = zeros(nofIt+1,1);
x_dense = zeros(nofBdryPts,2,(nofObs-1)*N);
alpha_dense = zeros(nofBdryPts,2,(nofObs-1)*N);

% obtaining initial estimate for theta based on the sparse observations
theta_hat(1,1) = theta_estimate_OU_2D(bdryPts_t,mu,dt_large);
theta_hat(1,1)

%theta_hat(:,:,1) = 2*ones(size(dummy(:,:,end)));

for it = 1:nofIt
  disp(sprintf('Iteration %i',it))
  total = 0;
  s = zeros(nofBdryPts,2,sampleSize);
  logL = zeros(sampleSize,1);
  
  
  nom = zeros(sampleSize,1);
  denom = zeros(sampleSize,1);
  
  for k = 1:sampleSize  
    % disp(sprintf('Sample %i',k))
    % simulating the proposal paths between the discrete observations
    for t = 1:nofObs-1 
        % sampling from the diffusion bridge
        x_dense(:,:,(t-1)*N+1) = bdryPts_t(:,:,t);
        
        % sampling from proposal distribution
        switch type
        case 'driftless'
	  for j=1:N-1 
	    K = K_matrix(x_dense(:,:,(t-1)*N+j),x_dense(:,:,(t-1)*N+j),sigma^2);
	    correction = (bdryPts_t(:,:,t) - x_dense(:,:,(t-1)*N+j))/(N*dt_small - j*dt_small);
	    x_dense(:,:,(t-1)*N+j+1) = x_dense(:,:,(t-1)*N+j) + dt_small*correction + sqrt(dt_small)*K*randn(size(mu));
	  end
	case 'drift'
	  for j=1:N-1  
	    K = K_matrix(x_dense(:,:,(t-1)*N+j),x_dense(:,:,(t-1)*N+j),sigma^2);
	    drift = theta_hat(it,1)*(mu - x_dense(:,:,(t-1)*N+j));
	    correction = (bdryPts_t(:,:,t) - x_dense(:,:,(t-1)*N+j))/(N*dt_small - j*dt_small);
	    x_dense(:,:,(t-1)*N+j+1) = x_dense(:,:,(t-1)*N+j) + dt_small*drift + dt_small*correction  + sqrt(dt_small)*K*randn(size(mu));
	  end
	end
	
	
	 
	% exact sampling
	%x_dense(:,(t-1)*N+1:t*N) = OU_bridge(bdryPts_t(:,t),bdryPts_t(:,t+1),mu,theta_hat(it,1),dt_small,dt_large,1);
        % likelihood(x_dense(:,(t-1)*N+1:t*N),theta_hat(it,1),mu,dt_small);
	logL(k,1) = logL(k,1) + likelihood(x_dense(:,:,(t-1)*N+1:t*N),theta_hat(it,1),mu,sigma,dt_small,type);
    end
  
%        	plot(x_dense)
%    	size(x_dense)
%    	pause()
    % L(k) = 1;
    % plot([1:dt:T]',x_dense(1,:))
    [dummy, nom(k,1), denom(k,1)] = theta_estimate_OU_2D(x_dense,mu,dt_small);
    % disp(dummy)
    
    % disp(norm(dummy(:,:,end) - theta_hat(:,:,1)))   
  end
  m = max(logL);
  L = exp(logL - m);
  
  theta_hat(it+1) = nom'*L/(denom'*L);
  %theta_hat(it+1) = sum(nom)/sum(denom);
  disp(sprintf('The estimate for theta is %d',theta_hat(it+1)))
  


  
end
end

function logL = likelihood(x_t,theta,mu,sigma,dt,type)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  logL = likelihood(x_t,theta,mu,dt,type)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


N = size(x_t,3);
dx_t = circshift(x_t,[0 0 -1]) - x_t;
T = N*dt;

% calculating the likelihood along the path
logL = 0;

for k = 2:N
  v = dx_t(:,:,k-1);
  x = x_t(:,:,k);
  % is this really the final value?
  x_final = x_t(:,:,end);
  drift = theta*(x_t(:,:,k) - mu);
  
  % calculate the drift
%    term1 = drift(:)'*v(:);
%    term2 = drift(:)'*drift(:);
%    term3 = (x(:) - x_final(:))'*(x(:)-x_final(:))/(2*(T - dt*(k-1)));
%    term4 = (x(:) - x_final(:))'*drift(:)/(T - dt*(k-1));
  % dont know if the third term is -/+ 1/2?
  
  K = K_matrix(x_t(:,:,k),x_t(:,:,k),sigma^2);
  K_inv = inv(K);
  K_inv_old = inv(K_matrix(x_t(:,:,k-1),x_t(:,:,k-1),sigma^2));
  
  term1 = drift(:)'*blkdiag(inv(K),inv(K))*v(:);
  term2 = drift(:)'*blkdiag(inv(K),inv(K))*drift(:);
  term3 = (x(:) - x_final(:))'*blkdiag(K_inv-K_inv_old,K_inv - K_inv_old)*(x(:)-x_final(:))/(2*(T - dt*(k-1)));
  term4 = (x(:) - x_final(:))'*(blkdiag(K_inv,K_inv))*drift(:)/(T - dt*(k-1));
  
  switch type
  case 'drift'
   % proposal with drift 
   logL = logL  - term3 - term4;
  case 'driftless'
    % proposal without drift
    logL = logL + term1 - term2*dt/2 - term3;
  end
   %logL
end

end

function [theta, nom, denom] = theta_estimate_OU_2D(x_t,mu,dt)
N = size(x_t,3);
T = round(N*dt);
dx_t = circshift(x_t,[0 0 -1]) - x_t;
%  disp(dx_t(:,1))
  
nom = 0;
denom = 0;

for k = 2:N
  v = dx_t(:,:,k-1);
  drift = (mu - x_t(:,:,k-1));
  
  % calculate the drift
  nom = nom + drift(:)'*v(:);
  denom = denom + drift(:)'*drift(:)*dt;
%    disp(nom)
%    disp(denom)

end
theta = nom/denom;
end

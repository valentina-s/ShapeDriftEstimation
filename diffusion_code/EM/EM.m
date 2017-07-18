function theta_hat = EM(bdryPts_t,sigma,dt,T,sampleSize,nofIt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% theta_hat = EM(bdryPts_t,sigma,dt,T,sampleSize,nofIt)
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nofBdryPts = size(bdryPts_t,1);
% the number of dense observations to generate between any two discrete observations
N = T/dt; 
nofObs = size(bdryPts_t,3);
theta_hat = zeros(nofBdryPts,2,nofIt+1);
x_dense = zeros(nofBdryPts,2,(nofObs-1)*N);
alpha_dense = zeros(nofBdryPts,2,(nofObs-1)*N);

% obtaining initial estimate for theta based on the sparse observations
dummy = drift_estimate(bdryPts_t,[],1,10);
theta_hat(:,:,1) = dummy(:,:,end);
%theta_hat(:,:,1) = 2*ones(size(dummy(:,:,end)));

for it = 1:nofIt 
  disp(sprintf('Iteration %i',it))
  total = 0;
  s = zeros(nofBdryPts,2,sampleSize);
  L = ones(sampleSize,	1);
  
  
  for k = 1:sampleSize  
    % disp(sprintf('Sample %i',k))
    % simulating the proposal paths between the discrete observations
    for t = 1:nofObs-1 
        [temp1 alpha_t] = Diffusion_bridge(bdryPts_t(:,:,t),bdryPts_t(:,:,t+1),theta_hat(:,:,it),sigma,dt,T);
	x_dense(:,:,(t-1)*N+1:t*N)  = temp1;
	alpha_dense(:,:,(t-1)*N+1:t*N) = alpha_t;
	L(k) = L(k) * pathLikelihood(temp1,[],theta_hat(:,:,it),sigma,T);
    end
    % L(k) = 1;
    
  
    [dummy s(:,:,k)] = drift_estimate(x_dense,[],dt,sigma);
    % disp(dummy)
    
    % disp(norm(dummy(:,:,end) - theta_hat(:,:,1)))   
  end
  % disp(L)
  temp = reshape(s,2*nofBdryPts,sampleSize)*L/(sum(L)*N);
  theta_hat(:,:,it+1) = reshape(temp,nofBdryPts,2);
  theta_hat(:,:,it+1)
  L
end
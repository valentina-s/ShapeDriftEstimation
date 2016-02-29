function L = pathLikelihood(x_t,alpha_t,theta,sigma,T)
N = size(x_t,3);
dt = T/N;

if isempty(alpha_t)
    dx_t = circshift(x_t,[0 0 -1]) - x_t;
else
    dx_t = alpha_t;
end

% calculating the likelihood along the path

logL = 0;
for k = 2:N
  v = dx_t(:,:,k-1);
  x = x_t(:,:,k);
  % is this really the final value?
  x_final = x_t(:,:,end);
  
  K = K_matrix(x_t(:,:,k),x_t(:,:,k),sigma^2);
  % K = eye(size(K));
  K_inv = inv(K);
  K_inv_old = inv(K_matrix(x_t(:,:,k-1),x_t(:,:,k-1),sigma^2));
  % K_inv_old = eye(size(K));
  drift = chol(K)*theta;
  
  % calculate the drift
  term1 = drift(:)'*blkdiag(inv(K),inv(K))*v(:);
  term2 = drift(:)'*blkdiag(inv(K),inv(K))*drift(:);
  term3 = (x(:) - x_final(:))'*blkdiag(K_inv-K_inv_old,K_inv - K_inv_old)*(x(:)-x_final(:))/(2*(T - dt*(k-1)));
  term4 = (x(:) - x_final(:))'*(blkdiag(K_inv,K_inv))*drift(:)/(T - dt*(k-1));
  % dont know if the third term is -/+ 1/2?
  
  % proposal with drift 
  logL = logL  - term3/2 - term4*dt;
  % proposal without drift
  % logL = logL + term1 - term2*dt/2 - term3;
end
logL
L = exp(logL-35)

function theta = theta_estimate_OU(x_t,mu,dt)
N = size(x_t,2);
T = round(N*dt);
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

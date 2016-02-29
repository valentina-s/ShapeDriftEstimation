function [m_new x_new alpha_new] = OU_update(bdryPts,ctrlPts,bdryPts_T,im,Sigma,Sigma_0,theta,dt,alpha)


nofCtrlPts = size(ctrlPts,1);

K_0 = K_matrix(ctrlPts,ctrlPts,Sigma_0^2);
K = K_matrix(ctrlPts,ctrlPts,Sigma^2 );
R = chol(K)/K_0;

alpha_noise = R'*randn(nofCtrlPts,2);

% calculate the gradient of the new state with respect to alpha
[G m_new x_new] = gradPhi(bdryPts,ctrlPts,Sigma_0,0.1,10,alpha);
BW_T = roipoly(im,bdryPts_T(:,1),bdryPts_T(:,2));

imagesc(BW_T);
par_obs = [1 0 1 1]';
E = - logObsL_region(BW_T,m_new,par_obs,'off');

% calculate the gradient of the distance with respect to the new state
g = gradU(m_new,bdryPts_T,im);
gg = reshape(G'*g(:),nofCtrlPts,2);

alpha_new =   -dt*theta*gg + sqrt(dt)*alpha_noise;

normGrad = alpha_new(:)'*blkdiag(K_0,K_0)*alpha_new(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Auxiliary functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function g = gradU(bdryPts,bdryPts_T,im)
BW = roipoly(im,bdryPts(:,1),bdryPts(:,2));
BW_T = roipoly(im,bdryPts_T(:,1),bdryPts_T(:,2));
% parameterizing by arc length, and calculating the normal
step = 1/size(bdryPts,1);
[m, nu] = arcLength(bdryPts, step);
f = computeF(BW_T,m);
g = grad(m,nu,f);

function g = grad(m,nu,f)
    F = interp2(f, m(:,1), m(:,2), 'linear', max(max(f)));
    g(:,1) = (F.*nu(:,1)); 
    g(:,2) = (F.*nu(:,2));  
end

function f = computeF(frame,bdryPts)  
  % setting likelihood parameters
  cin = 1; cout = 0;
  varin = 1; varout = 1;
  
  % Calculating the energy inside and outside

  E_1 = (frame - cin).^2/(2*varin);
  E_2 = (frame - cout).^2/(2*varout);

  f = E_1 - E_2;
end

end

end

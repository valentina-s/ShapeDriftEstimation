function [bdryPts_t ctrlPts_t alpha_t] = Diffusion_drift_OU(bdryPts,ctrlPts,bdryPts_T,im,sig_0,sig,theta,dt,T)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function [bdryPts_t ctrlPts_t alpha_t] = Diffusion_drift_OU(bdryPts,ctrlPts,bdryPts_T,im,sig_0,sig,theta,dt,T)
%
% Diffusion_drift_OU is a function which generates a diffusion with
% Ornstein-Uhlenbeck drift: i.e. a gradient of the distance to some fixed
% shape.
% theta \nabla dist(shape - mean_shape). 
%
% INPUTS:
% bdryPts
% ctrlPts
% bdryPts_T - template shape
% im - background image
% sig_0 - the deformation kernel width
% sig - the covariance kernel width
% theta
% dt
% T
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% this parameter determines whether to plot the contours
plotting = true;

N = T/dt;

nofBdryPts = size(bdryPts,1);
nofCtrlPts = size(ctrlPts,1);
bdryPts_t = zeros(nofBdryPts,2,N+1);
bdryPts_t(:,:,1) = bdryPts;
ctrlPts_t = zeros(nofCtrlPts,2,N+1);
ctrlPts_t(:,:,1) = ctrlPts;

alpha_t = zeros(nofCtrlPts,2,N);


% figure(1)
hold off
h1 = plot(bdryPts(:,1),bdryPts(:,2),'r--');
hold on
h2 = plot(ctrlPts(:,1),ctrlPts(:,2),'ro');
hold on

axis([0 100 0 100])
pause(0.1)


for i = 1:N
    
    % calculating the kernel matrices
    K_bdry = K_matrix(bdryPts_t(:,:,i),ctrlPts_t(:,:,i),sig_0^2);
    K_ctrl = K_matrix(ctrlPts_t(:,:,i),ctrlPts_t(:,:,i),sig_0^2);
   


    % calculate drift term with respect to the distance function
    g = gradU(bdryPts_t(:,:,i),bdryPts_T,im);
    % projecting the gradient
    g1  = K_bdry'*g;

    %max(max(g))
    %g
    % obtaining the corresponging momentum
    galpha = K_ctrl\g1;
    %max(max(galpha))

    % generate a random tangent vector at x
    C = K_matrix(ctrlPts_t(:,:,i),ctrlPts_t(:,:,i),sig*sig);
    K = K_matrix(ctrlPts_t(:,:,i),ctrlPts_t(:,:,i),sig_0*sig_0);
   
    R  = chol(C)/K;
    % noise = inv(chol(K))*alpha_drift + sqrt(dt)*inv(chol(K))*randn(size(alpha_drift));
    alpha_t(:,:,i) = - dt*theta*galpha  + sqrt(dt)*R'*randn(nofCtrlPts,2);
    % move along the exponential map
    bdryPts_t(:,:,i+1) = exponentialR(bdryPts_t(:,:,i),ctrlPts_t(:,:,i),alpha_t(:,:,i),sig_0,0.1,10);
    ctrlPts_t(:,:,i+1) = exponentialR(ctrlPts_t(:,:,i),ctrlPts_t(:,:,i),alpha_t(:,:,i),sig_0,0.1,10);


    if plotting
    
        set(h1,'Color','b')
        set(h2,'Color','b')

        h1 = plot(bdryPts_t(:,1,i+1),bdryPts_t(:,2,i),'r','Linewidth',2);
        h2 = plot(ctrlPts_t(:,1,i+1),ctrlPts_t(:,2,i),'ro','Linewidth',2);
   
        pause(0.01)
    end
end
end

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

function g = grad(m,nu,f)
    F = interp2(f, m(:,1), m(:,2), 'linear', max(max(f)));
    g(:,1) = (F.*nu(:,1)); 
    g(:,2) = (F.*nu(:,2));  
end
    

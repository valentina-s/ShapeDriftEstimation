function [bdryPts_t ctrlPts_t alpha_t] = Diffusion(bdryPts,ctrlPts,alpha_drift,sigma,dt,T)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% [bdryPts_t ctrlPts_t alpha_t] = Diffusion(bdryPts,ctrlPts,alpha_drift,sigma,dt,T)
%
% x0 - the initial configuration of the landmarks
% drift - contains the coefficients of the
% sigma - the deformation kernel width
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% this parameter determines whether to plot the contours
plotting = true;

N = T/dt;

nofCtrlPts = size(ctrlPts,1);
nofBdryPts = size(bdryPts,1);
bdryPts_t = zeros(nofBdryPts,2,N+1);
ctrlPts_t = zeros(nofCtrlPts,2,N+1);
bdryPts_t(:,:,1) = bdryPts;
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

    % generate a random tangent vector at x
    K = K_matrix(ctrlPts_t(:,:,i),ctrlPts_t(:,:,i),sigma*sigma);
    % K = eye(size(K));
    alpha_t(:,:,i) = dt*alpha_drift + sqrt(dt)*inv(sqrtm(K))*randn(size(alpha_drift));
    % move along the exponential map
    bdryPts_t(:,:,i+1) = exponentialR(bdryPts_t(:,:,i),ctrlPts_t(:,:,i),alpha_t(:,:,i),sigma,1,1);
    ctrlPts_t(:,:,i+1) = exponentialR(ctrlPts_t(:,:,i),ctrlPts_t(:,:,i),alpha_t(:,:,i),sigma,1,1);


    if plotting

        set(h1,'Color','b')
        set(h2,'Color','b')

        h1 = plot(bdryPts_t(:,1,i+1),bdryPts_t(:,2,i),'r--','Linewidth',2);
        h2 = plot(ctrlPts_t(:,1,i+1),ctrlPts_t(:,2,i),'ro','Linewidth',2);
   
        pause(0.01)

    end
end

    
    
    

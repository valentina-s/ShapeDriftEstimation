function [bdryPts_t ctrlPts_t alpha_t] = Diffusion_drift_dist(bdryPts,ctrlPts,bdryPts_T,sig_0,sig,theta,dt,T)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function [bdryPts_t ctrlPts_t alpha_t] = Diffusion_drift_LA(bdryPts,ctrlPts,meanArea,meanLength,delta_1,delta_2,sig_0,sig,dt,T)
%
% Diffusion_drift_LA is a function which evolves a shape according to a Langevin-type diffusion process.
% In this version the potential function U(shape) consist of two terms:
% U_1 = \delta_1*(length(shape) - length(template))^2
% U_2 = \delta_2*(area(shape) - area(template))^2,
% and its gradient is included as a drift term. 
%
%
% x0 - the initial configuration of the landmarks
% drift - contains the coefficients of the
% sig_0 - the deformation kernel width
% sig - the covariance kernel width
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% this parameter determines whether to plot the contours
plotting = false;

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
    g =  bdryPts_T - bdryPts_t(:,:,i);
    % projecting the gradient
    g1  = K_bdry'*g;

    % get the gradient directly with respect to the control points
    ctrlP

    % obtaining the corresponging momentum
    galpha = K_ctrl\g1;

    % generate a random tangent vector at x
    C = K_matrix(ctrlPts_t(:,:,i),ctrlPts_t(:,:,i),sig*sig);
    K = K_matrix(ctrlPts_t(:,:,i),ctrlPts_t(:,:,i),sig_0*sig_0);
   
    R  = chol(C)/K;
    % noise = inv(chol(K))*alpha_drift + sqrt(dt)*inv(chol(K))*randn(size(alpha_drift));
    alpha_t(:,:,i) =  dt*theta*galpha  + sqrt(dt)*R'*randn(nofCtrlPts,2);
    % move along the exponential map
    bdryPts_t(:,:,i+1) = exponentialR(bdryPts_t(:,:,i),ctrlPts_t(:,:,i),alpha_t(:,:,i),sig_0,1,1);
    ctrlPts_t(:,:,i+1) = exponentialR(ctrlPts_t(:,:,i),ctrlPts_t(:,:,i),alpha_t(:,:,i),sig_0,1,1);

    


    if plotting
    
        set(h1,'Color','b')
        set(h2,'Color','b')

        h1 = plot(bdryPts_t(:,1,i+1),bdryPts_t(:,2,i),'r--','Linewidth',2);
        h2 = plot(ctrlPts_t(:,1,i+1),ctrlPts_t(:,2,i),'ro','Linewidth',2);
   
        pause(0.01)
    end
end
end





    

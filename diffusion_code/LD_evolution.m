function [m x alpha] = LD_evolution(bdryPts,ctrlPts,meanArea,meanLength,image,Sigma,Sigma_0,dt,nofSteps,delta_1,delta_2,nofIt) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function [m x alpha] = LD_evolution(bdryPts,ctrlPts,meanArea,meanLength,image,Sigma,Sigma_0,dt,nofSteps,delta_1,delta_2,nofIt) 
%
% LD_evolution is a function which evolves a shape according to a Langevin diffusion process.
% In this version the potential function U(shape) consist of two terms:
% U_1 = \delta_1*(length(shape) - length(template))^2
% U_2 = \delta_2*(area(shape) - area(template))^2,
% and its gradient is included as a drift term. 
% 
% Uses LD_update.m
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Setting the dimensions
nofBdryPts = size(bdryPts,1);
nofCtrlPts = size(ctrlPts,1);

% Initializing some variables
m = zeros(nofBdryPts,2,nofSteps);
x = zeros(nofCtrlPts,2,nofSteps);
alpha = zeros(nofCtrlPts,2,nofSteps);

m(:,:,1) = bdryPts;
x(:,:,1) = ctrlPts;
alpha(:,:,1) = zeros(nofCtrlPts,2);


figure(1)

% OU evolution
for it = 1:nofIt
  
  if (it>1)
    
    [m(:,:,it) x(:,:,it)] = reparam(m(:,:,it),x(:,:,it));
    % project alpha on new control points
    alpha(:,:,it) = K_matrix(x(:,:,it),x(:,:,it),Sigma_0)\(K_matrix(x(:,:,it),x(:,:,it-1),Sigma_0)*alpha(:,:,it));
  end

  % update bdryPts, ctrlPts, and alpha
  [m(:,:,it+1) x(:,:,it+1) alpha(:,:,it+1)] = LD_update(m(:,:,it),x(:,:,it),meanArea,meanLength,image,Sigma,Sigma_0,dt,nofSteps,delta_1,delta_2,nofIt,alpha(:,:,it));


  % plotting
  imagesc(image)
  hold on
  plot(m(:,1,it+1),m(:,2,it+1),'.','Linewidth',2)
  hold off
  pause(0.1)

end
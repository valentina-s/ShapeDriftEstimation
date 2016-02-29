function [m x alpha] = OU_evolution(bdryPts,ctrlPts,bdryPts_T,im,Sigma,Sigma_0,theta,dt,T) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% [m x alpha] = OU_evolution(bdryPts,ctrlPts,bdryPts_T,ctrlPts_T,image,Sigma,Sigma_0,dt,nofSteps,delta,nofIt) 
%
% OU_evolution is a function which evolves a shape according the OU_evolutoin process on shape space.
% The distance to the mean is represented by a mismatch function between the two curves. 
% 
% Uses OU_update.m
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% determining whether to display the graphs
plotting = false;

% Determining the number of steps

nofSteps = T/dt;

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
hold off
% OU evolution
for it = 1:nofSteps

    [m(:,:,it+1) x(:,:,it+1) alpha(:,:,it+1)] = OU_update(m(:,:,it),x(:,:,it),bdryPts_T,im,Sigma,Sigma_0,theta,dt,alpha(:,:,it));


    if plotting

        hold on
        plot(bdryPts_T(:,1),bdryPts_T(:,2),'r','Linewidth',2)
        plot(m(:,1,it+1),m(:,2,it+1),'Linewidth',2)
        plot(x(:,1,it+1),x(:,2,it+1),'ro','Linewidth',2)
        hold off
        pause(0.1)
   end 
end

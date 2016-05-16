function [bdryPts_t ctrlPts_t alpha_t] = Diffusion_drift_LA(bdryPts,ctrlPts,meanArea,meanLength,theta,sig_0,sig,dt,T)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function [bdryPts_t ctrlPts_t alpha_t] = Diffusion_drift_LA(bdryPts,ctrlPts,meanArea,meanLength,theta,sig_0,sig,dt,T)
%
% Diffusion_drift_LA is a function which evolves a shape according to a Langevin-type diffusion process.
% In this version the potential function U(shape) consist of two terms:
% U_1 = \delta_1*(length(shape) - length(template))^2
% U_2 = \delta_2*(area(shape) - area(template))^2,
% and its gradient is included as a drift term. 
%
% theta = [area coefficient length coefficient]
%
% x0 - the initial configuration of the landmarks
% drift - contains the coefficients of the
% sig_0 - the deformation kernel width
% sig - the covariance kernel width
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% extracting the coefficients
delta_1 = theta(1); % area coefficient 
delta_2 = theta(2); % length coefficient

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
   

    % calculate drift term with respect to area
    [g_1 A] = gradAreaTerm(bdryPts_t(:,:,i),ctrlPts_t(:,:,i));
    % projecting the gradient
    g_1 = K_bdry'*g_1;
    
    galpha_1(:,1) = K_ctrl\g_1(:,1);
    galpha_1(:,2) = K_ctrl\g_1(:,2);

    % calculate drift term with respect to length
    % returns directly the projected gradient
    [g_2 L] = gradLengthTerm(bdryPts_t(:,:,i),ctrlPts_t(:,:,i),sig_0);
    
    galpha_2(:,1) = K_ctrl\g_2(:,1);
    galpha_2(:,2) = K_ctrl\g_2(:,2);
    

    % disp(sprintf('Area = %d; Length = %d',A,L))

    % generate a random tangent vector at x
    C = K_matrix(ctrlPts_t(:,:,i),ctrlPts_t(:,:,i),sig*sig);
    K = K_matrix(ctrlPts_t(:,:,i),ctrlPts_t(:,:,i),sig_0*sig_0);
   
    R  = chol(C)/K;
    % noise = inv(chol(K))*alpha_drift + sqrt(dt)*inv(chol(K))*randn(size(alpha_drift));
    alpha_t(:,:,i) = - dt*(delta_1*(A-meanArea)*galpha_1 + delta_2*(L - meanLength)*galpha_2)  + sqrt(dt)*R'*randn(nofCtrlPts,2);
    % move along the exponential map
    bdryPts_t(:,:,i+1) = exponentialR(bdryPts_t(:,:,i),ctrlPts_t(:,:,i),alpha_t(:,:,i),sig_0,0.1,10);
    ctrlPts_t(:,:,i+1) = exponentialR(ctrlPts_t(:,:,i),ctrlPts_t(:,:,i),alpha_t(:,:,i),sig_0,0.1,10);


    if plotting
    
        set(h1,'Color','b')
        set(h2,'Color','b')

        h1 = plot(bdryPts_t(:,1,i+1),bdryPts_t(:,2,i),'r--','Linewidth',2);
        h2 = plot(ctrlPts_t(:,1,i+1),ctrlPts_t(:,2,i),'ro','Linewidth',2);
   
        pause(0.01)
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Auxiliary functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [g length] = gradLengthTerm(bdryPts, ctrlPts, Sigma_0)
  nofBdryPts = size(bdryPts,1);
  nofCtrlPts = size(ctrlPts,1);


  Tangents = (circshift(bdryPts,-1) -circshift(bdryPts,1))/2;
  normTangents = sqrt(sum(Tangents.*Tangents,2));
  length = sum(normTangents);

  K_bdry = K_matrix(bdryPts,ctrlPts,Sigma_0*Sigma_0);

  for mm=1:nofCtrlPts
    s = zeros(1,2);
    for i=1:nofBdryPts
      dK_0 =  - (bdryPts(i,:) - ctrlPts(mm,:))*K_bdry(i,mm)/(Sigma_0*Sigma_0); 
      s = s + (Tangents(i,:)*dK_0')*Tangents(i,:)/normTangents(i); 
    end
    g(mm,:) = s;
  end 
end


function [g A] = gradAreaTerm(bdryPts, ctrlPts, Sigma_0)
    nofBdryPts = size(bdryPts,1);
    nofCtrlPts = size(ctrlPts,1);
    
    
    Tangents = (circshift(bdryPts,-1) -circshift(bdryPts,1))/2;
    g = [-Tangents(:,2) Tangents(:,1)];
    g = -g;
    BW = roipoly(zeros(100),bdryPts(:,1),bdryPts(:,2));
    A = sum(sum(BW));
end

function g = gradU(bdryPts,bdryPts_T,image)
  BW = roipoly(image,bdryPts(:,1),bdryPts(:,2));
  BW_T = roipoly(image,bdryPts_T(:,1),bdryPts_T(:,2));
  % parameterizing by arc length, and calculating the normal
  step = 1/size(bdryPts,1);
  [m, nu] = arcLength(bdryPts, step);
  f = computeF(BW_T,m);
  g = grad(m,nu,f);

end  
    

function theta_hat = drift_LA_estimate(x_t,xx_t,meanArea,meanLength,alpha_t,dt,sigma)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function theta_hat = drift_LA_estimate(x_t,xx_t,alpha_t,dt,sigma)
%
% x_t is the sequence of observations 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% rough version - dX = X_{t+1} - X_{t}
if isempty(alpha_t)
    dx_t = circshift(x_t,[0 0 -1]) - x_t;
else
    dx_t = alpha_t;
end
N = size(dx_t,3);


% calculate the average

% store the matrix in M and the vector in b
theta_hat = zeros(2,1);
M = zeros(2);
b = zeros(2,1);

theta_hat = zeros(2,N-1);

for i = 1:N-1


    % calculating the kernel matrices
    K_bdry = K_matrix(x_t(:,:,i),xx_t(:,:,i),sigma^2);
    K_ctrl = K_matrix(xx_t(:,:,i),xx_t(:,:,i),sigma^2);

    % calculate the gradient of the drift term with respect to area
    [gA A] = gradAreaTerm(x_t(:,:,i),xx_t(:,:,i));

    gA  = K_bdry'*gA; 

    galphaA(:,1) = K_ctrl\gA(:,1);
    galphaA(:,2) = K_ctrl\gA(:,2);

    % calculate drift term with respect to length
    % returns directly the projected gradient
    [gL L] = gradLengthTerm(x_t(:,:,i),xx_t(:,:,i),sigma);
    
    galphaL(:,1) = K_ctrl\gL(:,1);
    galphaL(:,2) = K_ctrl\gL(:,2);
    galphaL = gL;

    drift = [(L - meanLength)*gL(:)';(A - meanArea)*gA(:)'];
   

    
    
    M = M + drift*blkdiag(inv(K_ctrl),inv(K_ctrl))*drift'*dt;

    if isempty(alpha_t)
       alpha = K_bdry'\dx_t(:,:,i);
       temp = K_ctrl*alpha;
    else
       temp = K_ctrl*alpha_t(:,:,i);

    end


    b = b + drift*blkdiag(inv(K_ctrl),inv(K_ctrl))*temp(:);


    %if  isempty(alpha_t)
    %    K = chol(K_matrix(xx_t(:,:,i),xx_t(:,:,i),sigma^2));
    %    
    %     + K'\(x_t(:,:,i+1) - x_t(:,:,i));
    %else
    %    K = chol(K_matrix(xx_t(:,:,i),xx_t(:,:,i),sigma^2));
    % end
    theta_hat(:,i) = - M\b;
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

function [g A] = gradAreaTerm(bdryPts, ctrlPts)
    nofBdryPts = size(bdryPts,1);       
    Tangents = (circshift(bdryPts,-1) -circshift(bdryPts,1))/2;
    g = - [-Tangents(:,2) Tangents(:,1)];
    BW = roipoly(zeros(100),bdryPts(:,1),bdryPts(:,2));
    A = sum(sum(BW));
end




% equivalently
% theta_hat = (x_t(:,:,end) - x_t(:,:,1))/(dt*(N-1));

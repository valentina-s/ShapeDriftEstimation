function theta_hat = drift_OU_estimate(x_t,xx_t,x_mean,alpha_t,dt,sigma)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function theta_hat = drift_OU_estimate(x_t,xx_t,x_mean,alpha_t,dt,sigma)
%
%
% x_t is the sequence of observations 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initializing the background image
im = zeros(100);

% rough version - dX = X_{t+1} - X_{t}
if isempty(alpha_t)
    dx_t = circshift(x_t,[0 0 -1]) - x_t;
else
    dx_t = alpha_t;
end
N = size(dx_t,3);
theta_hat = zeros(N-1,1);


% calculate the average

% store the matrix in M and the vector in b

M = 0;
b = 0;
for i = 1:N-1


    % calculating the kernel matrices
    K_bdry = K_matrix(x_t(:,:,i),xx_t(:,:,i),sigma^2);
    K_ctrl = K_matrix(xx_t(:,:,i),xx_t(:,:,i),sigma^2);

    % calculate the gradient of the distance function
    g = gradU(x_t(:,:,i),x_mean,im);

    g1  = K_bdry'*g; 


    galpha(:,1) = K_ctrl\g1(:,1);
    galpha(:,2) = K_ctrl\g1(:,2);


    %cond(K_ctrl)

    drift = galpha(:);
    drift = g1(:);
    
    M = M + drift'*blkdiag(inv(K_ctrl),inv(K_ctrl))*drift*dt;

    if isempty(alpha_t)
       alpha = K_bdry'\dx_t(:,:,i);
       temp = K_ctrl*alpha;
    else
       temp = K_ctrl*alpha_t(:,:,i);
    end


    b = b + drift'*blkdiag(inv(K_ctrl),inv(K_ctrl))*temp(:);

    %if  isempty(alpha_t)
    %    K = chol(K_matrix(xx_t(:,:,i),xx_t(:,:,i),sigma^2));
    %    
    %     + K'\(x_t(:,:,i+1) - x_t(:,:,i));
    %else
    %    K = chol(K_matrix(xx_t(:,:,i),xx_t(:,:,i),sigma^2));
    % end
    theta_hat(i) = -b/M;
end

%theta_hat = - b/M;

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

% equivalently
% theta_hat = (x_t(:,:,end) - x_t(:,:,1))/(dt*(N-1));
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



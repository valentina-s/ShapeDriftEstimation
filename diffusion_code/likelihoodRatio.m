function L = likelihoodRatio(bdryPts_t,ctrlPts_t,alpha_t,theta,sigma,T,drift_type,par)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% likelihoodRatio calculates the likelihood ratio of the laws of two diffusions
% 
% INPUT:
%	bdryPts_t
%	ctrlPts_t
%	alpha_t
%	theta -	can be nx2-dimensional, 1-dimensional, 2-dimensional
%		(depending on the drift type), if we compare two processes, it
%		corresponds to the difference of the two parameters 
%		since the drift is linear wrt theta, we combine the two parameters in
%		in one difference
%	sigma -	the width of the kernel
%	T -	final time
%	drift_type - 'constant','OU','LA'
%	par - 	extra parameters for each drift
%      		par.bdryPts_mean - the template shape for OU drift
%		par.meanArea - the mean area for LA drift
%		par.meanLength - the mean length for LA drift
%		(constant drift does not have extra parameters)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





% this function calculates likelihood ratio based on two parameters
N = size(bdryPts_t,3);
dt = T/N;

if isempty(alpha_t)
    dx_t = circshift(ctrlPts_t,[0 0 -1]) - ctrlPts_t;
else
    dx_t = alpha_t;
end

% calculating the likelihood along the path
logL = 0;

for k = 1:N-1
  
  % constructing the kernels
  K_bdry = K_matrix(bdryPts_t(:,:,k),ctrlPts_t(:,:,k),sigma^2);
  K_ctrl = K_matrix(ctrlPts_t(:,:,k),ctrlPts_t(:,:,k),sigma^2);
  K_inv = inv(K_ctrl);




  switch drift_type
  
  
    case 'constant'
      drift = theta;

    case 'OU'

      im = zeros(100);
      % calculate the gradient of the distance function
      g = gradU(bdryPts_t(:,:,k),par.bdryPts_mean,im);
      g1  = K_bdry'*g; 


      galpha(:,1) = K_ctrl\g1(:,1);
      galpha(:,2) = K_ctrl\g1(:,2);

      drift = theta*galpha(:);
      % drift = g1(:);
  
    case 'LA'
      % calculate the gradient of the drift term with respect to area
      [gA A] = gradAreaTerm(bdryPts_t(:,:,k),bdryPts_t(:,:,k));

      gA  = K_bdry'*gA; 

      galphaA(:,1) = K_ctrl\gA(:,1);
      galphaA(:,2) = K_ctrl\gA(:,2);

      % calculate drift term with respect to length
      % returns directly the projected gradient
      [gL L] = gradLengthTerm(bdryPts_t(:,:,k),bdryPts_t(:,:,k),sigma);
      gL  = K_bdry'*gL; 

      K_ctrl\gL(:,1);
      galphaL(:,1) = K_ctrl\gL(:,1);
      galphaL(:,2) = K_ctrl\gL(:,2);
      %galphaL = gL;

      drift = theta(1,1)*(L - par.meanLength)*gL(:)'+ theta(2,1)*(A - par.meanArea)*gA(:)';
  end    


  v = dx_t(:,:,k);

  % calculate the drift
  term1 = drift(:)'*blkdiag(K_ctrl,K_ctrl)*v(:);
  term2 = drift(:)'*blkdiag(K_ctrl,K_ctrl)*drift(:);
  
  % proposal without drift
  logL = logL + term1 - term2*dt/2;
end

logL
L = exp(logL);
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


function [m_new x_new alpha_new] = LD_update(bdryPts,ctrlPts,mean_area,mean_length,image,Sigma,Sigma_0,dt,nofSteps,delta_1,delta_2,nofIt,alpha)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function [m_new x_new alpha_new] = LD_update(bdryPts,ctrlPts,mean_area,mean_length,image,Sigma,Sigma_0,dt,nofSteps,delta_1,delta_2,nofIt,alpha)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nofCtrlPts = size(ctrlPts,1);

% calculate the noise term
K_0 = K_matrix(ctrlPts,ctrlPts,Sigma_0^2);
K = K_matrix(ctrlPts,ctrlPts,Sigma^2 );
R = chol(K)*inv(K_0);

alpha_noise = R'*randn(nofCtrlPts,2);

% calculate drift term with respect to area
K_bdry = K_matrix(bdryPts,ctrlPts,Sigma_0^2);
[g area] = gradAreaTerm(bdryPts,ctrlPts,Sigma_0);
galpha_1 =  K_bdry\g;

% calculate drift term with respect to length
[g length] = gradLengthTerm(bdryPts,ctrlPts,Sigma_0);
galpha_2 =  K_0\g;


disp(sprintf('Area = %d; Length = %d',area,length))

% calculate the gradient of the new state with respect to alpha
[G m_new x_new] = gradPhi(bdryPts,ctrlPts,Sigma_0,dt,nofSteps,alpha);

% updating alpha
alpha_new =   1*alpha - 2*delta_1*(area-mean_area)*galpha_1 - 2*delta_2*(length - mean_length)*galpha_2 + 0*alpha_noise;


% normGrad = alpha_new(:)'*blkdiag(K_0,K_0)*alpha_new(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Auxiliary functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [g length] = gradLengthTerm(bdryPts, ctrlPts, Sigma_0)
  nofBdryPts = size(bdryPts,1);
  nofCtrlPts = size(ctrlPts,1);


  Tangents = (circshift(bdryPts,-1) -circshift(bdryPts,1))/2;
  normTangents = sqrt(sum(Tangents.*Tangents,2));
  length = sum(normTangents);

  K_bdry = K_matrix(bdryPts,ctrlPts,Sigma_0);

  for mm=1:nofCtrlPts
    s = zeros(1,2);
    for i=1:nofBdryPts
      dK_0 =  - (bdryPts(i,:) - ctrlPts(mm,:))*K_bdry(i,mm)/(Sigma_0*Sigma_0); 
      s = s + (Tangents(i,:)*dK_0')*Tangents(i,:)/normTangents(i); 
    end
    g(mm,:) = s;
  end 
end


function [g area] = gradAreaTerm(bdryPts, ctrlPts, Sigma_0)
  nofBdryPts = size(bdryPts,1);
  nofCtrlPts = size(ctrlPts,1);


  Tangents = (circshift(bdryPts,-1) -circshift(bdryPts,1))/2;
  g = [-Tangents(:,2) Tangents(:,1)];
  g = -g;
  BW = roipoly(image,bdryPts(:,1),bdryPts(:,2));
  area = sum(sum(BW));
  % K_bdry = K_matrix(bdryPts,ctrlPts,Sigma_0);
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

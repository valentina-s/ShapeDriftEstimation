function [m_new x_i_new alpha_new]= exponentialR(m, x_i, alpha, Sigma_0, dt, nofSteps)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% [m_new x_i_new alpha_new]= exponentialR(m, x_i, alpha, Sigma_0, dt, nofSteps)
% 
% This function evolves a set of boundary points m and a set of cotrol points x_i
% along the sub-Riemannian exponential map with initial momentum alpha. The width of the kernel is Sigma_0,
% dt is the size of the time step and nofSteps is the number of steps in the evolution.
% 
% 
% OUTPUT:
%       m_new - the final position of the curve
%       x_i_new - the final position of the control points
%	alpha_new - the new position of the momentum
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nofBdryPts = size(m,1);
nofCtrlPts = size(x_i,1);

%creating output variables
alpha_new = alpha;
  m_new = m;
  x_i_new = x_i;
for step = 1:nofSteps

  [dx dalpha] = updateODE(m,x_i_new,alpha_new,Sigma_0);

  % this is Laurent's version
  % [dx dalpha] = updateODE1(x_i_new,alpha_new,alpha_new,Sigma_0,size(alpha,1),2);
  
  % Updating the control points
  x_i_new = x_i_new + dt*dx;
  % Updating the boundary points
  m_new = m_new + dt*dm;
  % Updating the momenta
  alpha_new= alpha_new + dt*dalpha;
end

function [dx dalpha] = updateODE(m,x_i,alpha,Sigma_0);

    %Create the K_0 kernel matrix which moves the control points
    dist = bsxfun(@minus,x_i(:,1),x_i(:,1).').^2 + ...
           bsxfun(@minus,x_i(:,2),x_i(:,2).').^2;
    K_0_ctrl = exp(-dist/(2*Sigma_0*Sigma_0));

    %Create the k_0 kernel matrix which moves the boundary points
    dist = bsxfun(@minus,x_i(:,1),m(:,1).').^2 + ...
           bsxfun(@minus,x_i(:,2),m(:,2).').^2;
    K_0_bdry = exp(-dist/(2*Sigma_0*Sigma_0));

    % I need to construct a matrix of differences
    z2 = x_i(:,1)*ones(1,nofCtrlPts);
    m2 = ones(nofCtrlPts,1)*(x_i(:,1))';
    dist1 =  z2 - m2 ;

    z2 = x_i(:,2)*ones(1,nofCtrlPts);
    m2 = ones(nofCtrlPts,1)*(x_i(:,2))';
    dist2 =  z2 - m2 ;


    dist11 = dist1'.*K_0_ctrl'/(Sigma_0*Sigma_0)*alpha(:,1);
    dist22 = dist2'.*K_0_ctrl'/(Sigma_0*Sigma_0)*alpha(:,2);
    dist12 = dist1'.*K_0_ctrl'/(Sigma_0*Sigma_0)*alpha(:,2);
    dist21 = dist2'.*K_0_ctrl'/(Sigma_0*Sigma_0)*alpha(:,1);

    Dv = zeros(2,2,nofCtrlPts);
    Dv(1,1,:) = dist11;
    Dv(2,2,:) = dist22;
    Dv(1,2,:) = dist12;
    Dv(2,1,:) = dist21;


    % Update for the control points
    v = [];
    v(:,1) = K_0_ctrl'*alpha(:,1);
    v(:,2) = K_0_ctrl'*alpha(:,2);
    dx = v;

    %Update for the boundary points
    v = [];
    v(:,1) = K_0_bdry'*alpha(:,1);
    v(:,2) = K_0_bdry'*alpha(:,2);
    dm = v;

    %Update or the alphas
    dalpha(:,1) = -(squeeze(Dv(1,1,:)).*alpha(:,1)+squeeze(Dv(1,2,:)).*alpha(:,2));
    dalpha(:,2) = -(squeeze(Dv(2,1,:)).*alpha(:,1)+squeeze(Dv(2,2,:)).*alpha(:,2));
end


function [ux, za, zb] = updateODE1(x, a, b,sigma,N,dim)
% Y are all landmark information stacked:  
% dim*N for the landmarks (x)
% dim*N for the first momentum (a)
% dim*N for the second momentum (b)
% total : M = 3 * dim * N



%[x, a, b] = vectorToData(Y, N, dim) ;

sx2 = sum(x.^2, 2) ;
dist = sx2*ones(1, N) - 2*x*(x') + ones(N, 1)*sx2' ;

%dist = dist /(2*sigma*sigma) ;
gamma0 = exp(-dist/(2*sigma*sigma)) ;
%gamma00 = gamma0 + 1e-5*exp(-dist/(.5*sigma*sigma)) ;
gamma1 = - gamma0 / (2*sigma*sigma) ;
ux = gamma0 * a;
vx = gamma0 * b;

ss = sum(x.*ux, 2)*ones(1,N) ;
ss2 = x*(ux)' ;
xux = ss  - ss2 - ss2' + ss' ;
ss = sum(x.*vx, 2)*ones(1,N) ;
ss2 = x*(vx)' ;
xvx = ss  - ss2 - ss2' + ss' ;

scpaa = a*a' ;
scpab = a*b';
scpab = scpab + scpab' ;
gaa = gamma1 .* scpaa ;
gab = gamma1 .* scpab ;
gu = gamma1 .* xux; 
gv = gamma1 .* xvx; 
zb = conjGrad(@(x,g) g*x , (gv*a - gu*b), 20, gamma0) ; 
%zb = gamma00\(gv*a - gu*b);
za = - 2*(repmat(sum(gaa, 2), 1, dim).* x - gaa*x); 
zb =  zb - (repmat(sum(gab, 2), 1, dim).* x - gab*x); 
end


end

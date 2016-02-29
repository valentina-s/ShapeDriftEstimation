function K = K_matrix(x1,x2,var)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function constructs a matrix whose size is the size of
% x1 times size of x2. Its ij entry is the G (x_i-x_j) where 
% G is the Gaussian kernel. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%      d1 = size(x1,1);
%      d2 = size(x2,1);
%      z1 = x1(:,1).^2*ones(1,d2);
%      z2 = x2(:,1).^2*ones(1,d1);
%      dist =  z1 - 2* x1(:,1)*x2(:,1)' + z2' ;
%      z1 = x1(:,2).^2*ones(1,d2);
%      z2 = x2(:,2).^2 *ones(1,d1) ;
%      dist =  dist + z1 - 2*x1(:,2)*x2(:,2)' + z2' ;


  dist = bsxfun(@minus,x1(:,1),x2(:,1)').^2 + ...
     bsxfun(@minus,x1(:,2),x2(:,2)').^2;
  
  % the segmentation has a normalization factor?
  % K = exp(-dist/(2*sigma*sigma))/(2*pi*sigma*sigma);
  K = exp(-dist/(2*var));
end
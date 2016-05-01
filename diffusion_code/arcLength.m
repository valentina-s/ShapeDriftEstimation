function [z, normal] = arcLength(z0, step) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arcLength            Rediscretization of a Closed Curve by a Fixed Step 
%
% Call:  [z, normal] = arcLength(z0, step) 
% Input:
% 	z0 - an N by 2 array, where z0(:,1) and z0(:,2) are the x and y coordinates of a CLOSED curve
% 	step - a double indicating the size of the discretization step ( default: step = 0.01)
% Output: 
% 	z - an N by 2 array, where z(:,1) and z(:,2) are the x and y coordinates for the rediscretized curve
%       normal - an N by 2 array of normals to the rediscretized curve


if (nargin == 1)
  step = 0.01 ;
end ;

if (sum(sum(abs(z0(1,:)-z0(end,:))))<1e-7)
    z0 = z0(1:end-1, :) ;
end

% computes the arc-length
dz = circshift(z0,-1) - z0 ;
ds = sqrt(dz(:,1).^2 + dz(:, 2).^2) ;
if (min(min(ds)) < 1e-7)
    disp('same points');
end;
s = zeros(size(z0,1)+1, 1) ;
s(1) = 0 ;
s(2:size(z0,1)+1) = cumsum(ds) ; 

% subsamples the arc-length
z = zeros(size(z0)) ;
L = sum(sum(ds)) ;
s = s/L ;
z(1, :) = z0(1, :) ;
oldk = 1 ;
l=1 ;
t = 0:step:1 ;
z2(:,1) = interp1(s, [z0(:,1); z0(1,1)], t(1:end-1)) ;
z2(:,2) = interp1(s,[z0(:,2); z0(1, 2)], t(1:end-1)) ;
z=z2 ;
dz = (circshift(z,-1) - circshift(z,1))/2 ;
normal = [-dz(:,2), dz(:,1)] ;
n = normal(:,1).^2 + normal(:,2).^2 ;
normal = normal ./ repmat(sqrt(1e-10 + n), 1, 2) ;
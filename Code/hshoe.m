% w=hshoe(x,y,z,x1,y1,z1,Gamma,N)
%
%function to find the downwash velocity (w_ij) induced on a point x_i,y_i,z_i 
%by N horseshoe vortices,each of strength Gamma and with end points 
%at x1(j,1),y1(j,1),z1(j,1) and x1(j,2),y1(j,2),z1(j,2)
%
% OUTPUTS
% w - N x 1 array of induced velocities (+ve downwards) (m/s)
% INPUTS
% y - scalar y-position to estimate downwash at (m)
% y1 - N x 2 array containing horseshoe vortex endpoints (m)
%      make sure that y1(:,1)<y1(:,2) to get the right direction for W
%      for y1(:,1)<y<y1(:,2) the answer should be +ve
% Gamma - N x 1 array containing the strength of each vortex (m^2/s)
% N - Scalar defining number of vortices

function w=hshoe(y,y1,Gamma)
    
    w=Gamma.*(1./(y-y1(:,1))+1./(y1(:,2)-y))/(4*pi);
    
end
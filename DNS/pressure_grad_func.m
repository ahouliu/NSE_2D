%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%  INPUT:  p(x,y)
%  OUTPUT: p_x,p_y
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [px,py] = pressure_grad_func(p)

%% physical space
N       = size(p,1);
N_half = N/2;
%% spectral space
k1 = [0:N_half-1,0, -N_half+1:-1]';
k2 = [0:N_half-1,0, -N_half+1:-1]';
p_hat = fft2(p);
p_hat(1,1) = 0;
p_hat(N_half+1,:) = 0;
p_hat(:,N_half+1) = 0;
%% evaluating pressure gradient
px = 1i*k1'.*p_hat;
py = 1i*k2 .*p_hat;
px = real(ifft2(px));
py = real(ifft2(py));
end
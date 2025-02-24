%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    -Laplacian(p) = 2(u_y v_x       - u_x v_y)       - (f1_x-f2_y)
%                  = 2(psi_xy psi_xy - psi_xx psi_yy) - (f1_x-f2_y)
%
%  INPUT:  w(x,y)
%  OUTPUT: p(x,y)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p] = pressure_func(psi_hat,t,w_force,f1xf2y_hat)

%% physical space
N       = size(psi_hat,1);
N_half = N/2;
% p       = zeros(N,N);
%% spectral space
k1 = [0:N_half-1,0, -N_half+1:-1]';
k2 = [0:N_half-1,0, -N_half+1:-1]';
K12 = bsxfun(@(x,y) x.^2+y.^2,k1',k2);
K12_1                     = 1./K12;
K12_1(1, 1)               = 0;
K12_1(1, N_half+1)        = 0;
K12_1(N_half+1, 1)        = 0;
K12_1(N_half+1, N_half+1) = 0;
% external force
% [X,Y]   = meshgrid((0:N-1)'*(Lx/N), (0:N-1)'*(Ly/N)); %(y,x)
% [F1,F2] = f12_func(X,Y,ep1,ep2);
% F1      = 1i*k1'.*fft2(F1);
% F2      = 1i*k2 .*fft2(F2);
% f1xf2y_hat                     = F1+F2;
% f1xf2y_hat(abs(f1xf2y_hat)<1e-10) = 0;
% f1xf2y_hat                     = 1i*imag(f1xf2y_hat);
% f1xf2y_hat                    = sparse(f1xf2y_hat);
%% evaluating pressure field
temp12 = k1'.*k2 .*psi_hat;
temp11 = k1'.*k1'.*psi_hat;
temp22 = k2 .*k2 .*psi_hat;
p = 2*( conv_func(temp12,temp12) - conv_func(temp11,temp22) ) - f1xf2y_hat*g_func(t,w_force);
p = p.*K12_1;
p = real(ifft2(p));

end

function [W] = conv_func(U,V)

N        = size(U,1);
M        = 1.5*N;
N_half   = 0.5*N;
M_half   = 0.5*M;
UU   = zeros(M,M);
VV   = zeros(M,M);

%% insert into M*M
UU(       1:N_half +1,  1:N_half +1)         =  U(1:N_half +1,       1:N_half +1);         %(+,+)
UU(       1:N_half +1,   end-N_half +2:end)  =  U(1:N_half +1,       end-N_half +2:end);   %(+,-)
UU( end-N_half +2:end,   1:N_half +1)        =  U(end-N_half +2:end, 1:N_half +1);         %(-,+)
UU( end-N_half +2:end,   end-N_half +2:end)  =  U(end-N_half +2:end, end-N_half +2:end);   %(-,-)
VV(      1:N_half +1,  1:N_half +1)          =  V(1:N_half +1,       1:N_half +1);         %(+,+)
VV(      1:N_half +1,  end-N_half +2:end)    =  V(1:N_half +1,       end-N_half +2:end);   %(+,-)
VV(end-N_half +2:end,  1:N_half +1)          =  V(end-N_half +2:end, 1:N_half +1);         %(-,+)
VV(end-N_half +2:end,  end-N_half +2:end)    =  V(end-N_half +2:end, end-N_half +2:end);   %(-,-)
%% spectral -> physical -> spectral
W   = real(ifft2(UU)).* real(ifft2(VV));
W   = 2.25*fft2(W); % M^2/N^2
%% only N*N are reserved
W(N_half +2:N+1,:) = [];
W(:,N_half +2:N+1) = [];% k1,k2 = -N_half +1,...,N_half 
W(N_half +1,:)     = 0;
W(:,N_half +1)     = 0; % W_k1 OR W_k2 = 0

end
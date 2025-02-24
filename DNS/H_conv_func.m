%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Evaluating the nonlinear terms. 3/2 rules
%  M =(3/2)N 
%  M^2/N^2=9/4 = 2.25
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [W] = H_conv_func(w_hat,psi_hat)
%% pre-stage
N        = size(w_hat,1);
M        = 1.5*N;
N_half   = 0.5*N;
M_half   = 0.5*M;
ww_hat   = zeros(M,M);
ppsi_hat = zeros(M,M);
%% insert into M*M
ww_hat(       1:N_half +1,  1:N_half +1)         =  w_hat(1:N_half +1,1:N_half +1);               %(+,+)
ww_hat(       1:N_half +1,   end-N_half +2:end)  =  w_hat(1:N_half +1,end-N_half +2:end);         %(+,-)
ww_hat( end-N_half +2:end,   1:N_half +1)        =  w_hat(end-N_half +2:end,1:N_half +1);         %(-,+)
ww_hat( end-N_half +2:end,   end-N_half +2:end)  =  w_hat(end-N_half +2:end,end-N_half +2:end);   %(-,-)
ppsi_hat(      1:N_half +1,  1:N_half +1)        =  psi_hat(1:N_half +1,1:N_half +1);             %(+,+)
ppsi_hat(      1:N_half +1,  end-N_half +2:end)  =  psi_hat(1:N_half +1,end-N_half +2:end);       %(+,-)
ppsi_hat(end-N_half +2:end,  1:N_half +1)        =  psi_hat(end-N_half +2:end,1:N_half +1);       %(-,+)
ppsi_hat(end-N_half +2:end,  end-N_half +2:end)  =  psi_hat(end-N_half +2:end,end-N_half +2:end); %(-,-)
%% spectral -> physical -> spectral
kk1 = [0:1:M_half,-M_half+1:1:-1]';
kk2 = [0:1:M_half,-M_half+1:1:-1]';
W   =   real(ifft2( 1i*kk2 .* ww_hat )).* real(ifft2( 1i*kk1'.* ppsi_hat ))...
      - real(ifft2( 1i*kk1'.* ww_hat )).* real(ifft2( 1i*kk2 .* ppsi_hat ));
W   = 2.25*fft2(W); % M^2/N^2
%% only N*N are reserved
W(N_half +2:N+1,:) = [];
W(:,N_half +2:N+1) = [];% k1,k2 = -N_half +1,...,N_half 
W(N_half +1,:)     = 0;
W(:,N_half +1)     = 0; % W_k1 OR W_k2 = 0
end
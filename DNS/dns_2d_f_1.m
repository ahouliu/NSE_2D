%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%  Fourier-Galerkin spectral method 2d N-S, [0,Lx]*[0,Ly],[t0, T_max]
%   STORING ALL temporal snapshots!!!
%  velocity-vorticity formulations ( vorticity-stream --> only stream, may unconditionally unstable!!( frefer to David Gottlieb, Steven A. Orszag - Numerical Analysis of Spectral Methods - Theory and Applications-)
%  N-- grid number in [0,Lx].
%  ETD-RK4 FROM Fourth-order time-stepping for stiff PDEs. A. Kassam & L. N. Trefethen
%
%  [X,Y]=mehsgrid(x,y)---(y,x), not (x,y) !!!
%  k1 =N/2 OR k2 = N/2 --> set to be 0.
%  in MATLAB, fft --> the Galerkin spectrum * N (N = number of grids), i.e. u_hat * N.
%
%  PARAMATERS:
%  pras={Lx, Ly, N, nu}. d--dt,t \in [e0, e1], (f1,f2) - strengths of the forces
%  w--the initial values of \omega
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [W_data,Ek,Dk,Ik,w,psi] = dns_2d_f_1(paras,d,e0,e1,f1,f2,Tf,w,psi)
tic;

Lx      = paras{1};
Ly      = paras{2};
LxLy    = Lx*Ly;
N       = paras{3};
N_half  = N/2;
p_DFT   = 1/(N^2); % from DFT in MATLAB to DFT in  reality
nu      = paras{4};
if isinf(Tf)
    w_force = 0;
else
    w_force = 2*pi/Tf;
end
dt      = d;
t0      = e0;
t_max   = e1;
Ninteration = round(t_max/dt)-round(t0/dt);
if abs(t_max-t0-dt*Ninteration)>1e-10
    error('What happend?');
end
clear('paras','d','e0','e1','t_max','Tf');
%% spectral space
k1 = [0:N/2-1,0, -N_half+1:-1]';
k2 = [0:N/2-1,0, -N_half+1:-1]';
K12 = bsxfun(@(x,y) x.^2+y.^2,k1',k2);
K12_1                     = 1./K12;
K12_1(1, 1)               = 0;
K12_1(1, N_half+1)        = 0;
K12_1(N_half+1, 1)        = 0;
K12_1(N_half+1, N_half+1) = 0;

w_hat                = fft2(w);
w_hat(1, 1)          = 0;
w_hat(N_half+1, :)   = 0;
w_hat(:, N_half+1)   = 0;
if ~exist('psi','var')||isempty(psi)
    psi_hat = w_hat.*K12_1;
else
    psi_hat = fft2(psi);
end
psi_hat(1, 1)        = 0;
psi_hat(N_half+1, :) = 0;
psi_hat(:, N_half+1) = 0;

% external force
f12_hat  = f12hat_func(N,f1,f2);
clear('Lx','Ly','f1','f2');
%% storing some data
if dt>=0.01
    sample_interval = dt;
else
    sample_interval = 0.01;
end
interval = round(sample_interval/dt);
LEN = round(Ninteration/interval) + 1;

W_data(N,N,LEN)   = 0;
Ek(LEN,1)         = 0;
Dk(LEN,1)         = 0;
Ik(LEN,1)         = 0;

W_data(:,:,1) = w;
Ek(1)         = 0.5*LxLy*sum(sum( K12.*abs(psi_hat).^2 ))*p_DFT ^2;
Dk(1)         = (LxLy*nu)* norm(w_hat,'fro')^2*p_DFT ^2;
Ik(1)         = LxLy* sum(sum( f12_hat.*conj(psi_hat) ))*p_DFT^2*g_func(t0,w_force);
%% pre-computing 
E_f  = exp(-1*dt*nu*K12);
E2_f = exp(-0.5*dt*nu*K12);
M = 32; % number of points for complex means
r = exp( 1i*pi* ((1:M)-1/2)/M );
[K1,K2,Lhr] = meshgrid(k1,k2,r);
Lhr = -1*dt*nu*(K1.^2+K2.^2) + Lhr;
clear('K1','K2','r','k1','k2');
Q_f      = dt*real(sum( (exp(0.5*Lhr)-1)./Lhr             ,3)/M);
ETDRK_f1 = dt*real(sum( (-4-Lhr+exp(Lhr).*(4-3*Lhr+Lhr.^2))./(Lhr.^3)  ,3)/M);
ETDRK_f2 = dt*real(sum( (2+Lhr+exp(Lhr).*(-2+Lhr))./(Lhr.^3)           ,3)/M);
ETDRK_f3 = dt*real(sum( (-4-3*Lhr-Lhr.^2+exp(Lhr).*(4-Lhr))./(Lhr.^3)  ,3)/M);
%% update flow fields
for ind = 1:Ninteration
    t_now = t0+(ind-1)*dt;
    Nw = H_conv_func(w_hat,psi_hat) + f12_hat*g_func(t_now,w_force);
    ETDRK_a = E2_f.*w_hat + Q_f.*Nw;
    ETDRK_a(1,1) = 0;
    ETDRK_a(N_half+1,:) = 0;
    ETDRK_a(:,N_half+1) = 0;
    psi_hat = ETDRK_a.*K12_1;
    Na = H_conv_func(ETDRK_a,psi_hat) + f12_hat*g_func(t_now+0.5*dt,w_force);
    ETDRK_b = E2_f.*w_hat + Q_f.*Na;
    ETDRK_b(1,1) = 0;
    ETDRK_b(N_half+1,:) = 0;
    ETDRK_b(:,N_half+1) = 0;
    psi_hat = ETDRK_b.*K12_1;
    Nb = H_conv_func(ETDRK_b,psi_hat) + f12_hat*g_func(t_now+0.5*dt,w_force);
    ETDRK_c = E2_f.*ETDRK_a + Q_f.*(2*Nb-Nw);
    ETDRK_c(1,1) = 0;
    ETDRK_c(N_half+1,:) = 0;
    ETDRK_c(:,N_half+1) = 0;
    psi_hat = ETDRK_c.*K12_1;
    Nc = H_conv_func(ETDRK_c,psi_hat) + f12_hat*g_func(t_now+dt,w_force);
    w_hat = E_f.*w_hat + ETDRK_f1.*Nw + 2*ETDRK_f2.*(Na+Nb) + ETDRK_f3.*Nc;
    w_hat(1,1) = 0;
    w_hat(N_half+1,:) = 0;
    w_hat(:,N_half+1) = 0;
    psi_hat = w_hat.*K12_1;
    
    E2 = 0.5*LxLy*sum(sum( K12.*abs(psi_hat).^2 ))*p_DFT ^2;
	if isnan(E2)
        error('Explode at t=%g (dt=%g).\n ',t0+dt*ind,dt);
    end
    if mod(ind,interval)==0
        jnd = ind/interval;
         W_data(:,:,1+jnd) = real(ifft2(w_hat));
         Ek(jnd+1)         =  E2;
         Dk(jnd+1)         = (LxLy*nu)*norm(w_hat,'fro')^2*p_DFT ^2;
         Ik(jnd+1)         = LxLy*sum(sum( f12_hat.*conj(psi_hat) ))*p_DFT ^2* g_func(t0+ind*dt,w_force);
    end
end
%% return values(?)
w   = real(ifft2(w_hat));
psi = real(ifft2(psi_hat));
toc
end
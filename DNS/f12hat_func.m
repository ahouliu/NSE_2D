%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%  Returning the external forces in (k_1,k_2)
%
%  WARNING!!!
%      in MATLAB, it is (k2,k1), NOT (k1,k2), while using the index (n1,n2)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [f12_hat] = f12hat_func(N,ep1,ep2)

k0 = 4;
% k0 = 6;
%warning(['wave number of the force is ',num2str(w0)]);

%% DFT型

% Lx =2*pi;Ly=2*pi;
% [X,Y]   = meshgrid((0:N-1)'*(Lx/N), (0:N-1)'*(Ly/N)); %(y,x)

% f1 =   (ep1/4).*(1-cos(k0*X)).* (1-cos(k0*Y));
% f2 =   (ep2/4).*(1-cos(k0*X)).* (1-cos(k0*Y));
% k1 = [0:N/2-1,0, -N/2+1:-1]';
% k2 = [0:N/2-1,0, -N/2+1:-1]';
% f12_hat = 1i*( -1*k2.*fft2(f1) + k1'.*fft2(f2) );

%% 指定型

f12_hat = sparse(N,N);
f12_hat(1,1+k0)     = -(ep2*k0/8)*1i; % ( k0,0)
f12_hat(1,end-k0+1) =  (ep2*k0/8)*1i; % (-k0,0)
f12_hat(1+k0,1)     =  (ep1*k0/8)*1i; % (0, k0)
f12_hat(end-k0+1,1) = -(ep1*k0/8)*1i; % (0,-k0)

f12_hat(k0+1,k0+1)         = -((ep1-ep2)*k0/16)*1i; % ( k0, k0)
f12_hat(end-k0+1,end-k0+1) =  ((ep1-ep2)*k0/16)*1i; % (-k0,-k0)
f12_hat(end-k0+1,k0+1)     =  ((ep1+ep2)*k0/16)*1i; % ( k0,-k0)
f12_hat(k0+1,end-k0+1)     = -((ep1+ep2)*k0/16)*1i; % (-k0, k0)
f12_hat=f12_hat*N^2;

%% 换个外力

% Lx = 2*pi;
% Ly = 2*pi;
% [X,Y] = meshgrid((0:N-1)'*(Lx/N), (0:N-1)'*(Lx/N)); %(y,x)
% f1 = ep1*exp(-(Y-Ly/2).^2).*cos(k0*Y);
% f2 = ep2*exp(-(X-Lx/2).^2).*cos(k0*X);
% k1 = [0:N/2-1,0, -N/2+1:-1]';
% k2 = [0:N/2-1,0, -N/2+1:-1]';
% f12_hat = 1i*( -1*k2.*fft2(f1) + k1'.*fft2(f2) );
% f12_hat(abs(f12_hat)<1e-10)=0;
% INDEX = find(abs(real(f12_hat))<1e-10);
% f12_hat(INDEX) = 1i*imag(f12_hat(INDEX));
% INDEX = find(abs(imag(f12_hat))<1e-10);
% f12_hat(INDEX) = real(f12_hat(INDEX));

% Lx = 2*pi;
% Ly = 2*pi;
% [X,Y] = meshgrid((0:N-1)'*(Lx/N), (0:N-1)'*(Lx/N)); %(y,x)
% f1 = ep1*exp(-((X-Lx/2).^2 + (Y-Ly/2).^2).^(2)).*cos(k0*Y);
% f2 = ep2*exp(-((X-Lx/2).^2 + (Y-Ly/2).^2).^(2)).*cos(k0*X);
% k1 = [0:N/2-1,0, -N/2+1:-1]';
% k2 = [0:N/2-1,0, -N/2+1:-1]';
% f12_hat = 1i*( -1*k2.*fft2(f1) + k1'.*fft2(f2) );
% f12_hat(abs(f12_hat)<1e-10)=0;
% INDEX = find(abs(real(f12_hat))<1e-10);
% f12_hat(INDEX) = 1i*imag(f12_hat(INDEX));
% INDEX = find(abs(imag(f12_hat))<1e-10);
% f12_hat(INDEX) = real(f12_hat(INDEX));

% k0 = 16;
% Lx = 2*pi;
% Ly = 2*pi;
% [X,Y] = meshgrid((0:N-1)'*(Lx/N), (0:N-1)'*(Lx/N)); %(y,x)
% RR = sqrt((X-Lx/2).^2 + (Y-Ly/2).^2 );
% L = 0.3*min(Lx/2,Ly/2);
% f12 = @(r) (abs(r)<L).*exp(1+1./(r.^2-1)).*cos(k0*r) + (abs(r)>=L).*0;
% f12 = 0.5*sqrt(ep1^2+ep2^2)*f12(RR);
% f12_hat = fft2(f12);
% f12_hat(1,1) = 0;
% f12_hat(N/2+1,:) = 0;
% f12_hat(:,N/2+1) = 0;
% f12_hat(abs(f12_hat)<1e-10)=0;
% INDEX = find(abs(real(f12_hat))<1e-10);
% f12_hat(INDEX) = 1i*imag(f12_hat(INDEX));
% INDEX = find(abs(imag(f12_hat))<1e-10);
% f12_hat(INDEX) = real(f12_hat(INDEX));

% k0 = 8;
% Lx = 2*pi;
% Ly = 2*pi;
% [X,~] = meshgrid((0:N-1)'*(Lx/N), (0:N-1)'*(Ly/N)); %(y,x)
% 
% OMEGA0 = 2*pi/Lx;
% f12 = @(x) (abs(x-0.5*Lx)<=0.5*Lx/k0).*(1-cos(k0*OMEGA0*X-(k0-1)*pi)) + (abs(x-0.5*Lx)>0.5*Lx/k0).*0;
% f12 = ep1*f12(X);
% f12_hat = fft2(f12);
% f12_hat(1,1) = 0;
% f12_hat(N/2+1,:) = 0;
% f12_hat(:,N/2+1) = 0;
% f12_hat(abs(f12_hat)<1e-10)=0;
% INDEX = find(abs(real(f12_hat))<1e-10);
% f12_hat(INDEX) = 1i*imag(f12_hat(INDEX));
% INDEX = find(abs(imag(f12_hat))<1e-10);
% f12_hat(INDEX) = real(f12_hat(INDEX));

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%  Returning the external forces in (k_1,k_2)
%
%  WARNING!!!
%      in MATLAB, it is (k2,k1), NOT (k1,k2), while using the index (n1,n2)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [f12_hat] = f12hat_func(N)

%% type-1
% k0 = 16;
% L = 2*pi;
% 
% k1  = [0:N/2-1,0, -N/2+1:-1]';
% k2  = (0:N/2-1)';
% K12 = bsxfun(@(x,y) sqrt(x.^2+y.^2),k1',k2);
% 
% KL                      = 7*N/(15*L);
% fk                      = @(k) exp(-(k-0.3*N/L).^4).*(1+cos(k0*k))/2;
% f12_hat                 = zeros(N/2,N);
% f12_hat(abs(K12/L)< KL) = fk(K12(abs(K12/L)< KL)/L);
% f12_hat = f12_hat.*exp(1i*rand(size(f12_hat)));
% f12_hat                 = Hsym(f12_hat,N);
% 
% f12_hat = real(ifft2(f12_hat));
% f12_hat = [f12_hat(:,N/2+1:end),f12_hat(:,1:N/2)];
% f12_hat = [f12_hat(N/2+1:end,:);f12_hat(1:N/2,:)];
% f12_hat = fft2(f12_hat/max(f12_hat(:)));
%% type-2
% k0    = 8;
% L    = 2*pi;
% [X,Y] = meshgrid((0:N-1)'*(L/N), (0:N-1)'*(L/N)); %(y,x)
% RR    = sqrt( (X-L/2).^2+(Y-L/2).^2 );
% k1  = [0:1:N/2-1,0, -N/2+1:-1]';
% k2  = [0:1:N/2-1,0, -N/2+1:-1]';
% K_max = 3*N/8;
% K_filter = bsxfun(@(x,y) exp(-(x/K_max).^16).*exp(-(y/K_max).^16) ,k1',k2);
% 
% KL  = L/6;
% fx = @(r) (abs(r-L/5)<KL).*exp(1+1./( ((r-L/5)/KL).^2-1 )).*exp(-(r-L/5).^2).*cos(k0*(r-L/5)) + (abs(r-L/5)>=KL).*0;
% f12_hat = zeros(size(RR));
% f12_hat(abs(RR-L/5)<KL) = fx(RR(abs(RR-L/5)<KL));
% f12_hat = fft2(f12_hat).*K_filter;
% f12_hat(1,1) = 0;
% f12_hat(:,N/2+1) = 0;
% f12_hat = f12_hat(1:N/2,:);
%% type-3
k0    = 8;
L    = 2*pi;
[X,Y] = meshgrid((0:N-1)'*(L/N), (0:N-1)'*(L/N)); %(y,x)
RR    = sqrt( (X-L/2).^2+(Y-L/2).^2 );
k1  = [0:1:N/2-1,0, -N/2+1:-1]';
k2  = [0:1:N/2-1,0, -N/2+1:-1]';
K_max = 3*N/8;
K_filter = bsxfun(@(x,y) exp(-(x/K_max).^16).*exp(-(y/K_max).^16) ,k1',k2);

fx = @(r) exp(-r.^2).*sin(k0*r) ;
f12_hat = fft2(fx(RR)).*K_filter;
f12_hat(1,1) = 0;
f12_hat(:,N/2+1) = 0;
f12_hat = f12_hat(1:N/2,:);
%% 调整
f12_hat(abs(f12_hat)<1e-10)=0;
INDEX = find(abs(real(f12_hat))<1e-10);
f12_hat(INDEX) = 1i*imag(f12_hat(INDEX));
INDEX = find(abs(imag(f12_hat))<1e-10);
f12_hat(INDEX) = real(f12_hat(INDEX));
f12_hat = Hsym(f12_hat(1:N/2,:),N);

A_fac = 1/max(abs(reshape(real(ifft2(f12_hat)),[],1)));
f12_hat = A_fac*f12_hat;
end

function fk = Hsym(fk,N)

fk(1,1) = 0;
fk(:,N/2+1)     = 0;
fk(1,N/2+2:end) = fliplr(fk(1,2:N/2));
fk = [fk;zeros(1,N);flipud(fk(2:N/2,1)),rot90(fk(2:N/2,2:end),2)];

end
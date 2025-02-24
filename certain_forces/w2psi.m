%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%  w(x,y) to \psi(x,y)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [psi,u,v] = w2psi(w)

if size(w,1)~=size(w,2)
    error('NOT square');
end

N = size(w,1);
k1 = [0:1:N/2-1,0, -N/2+1:-1]';
k2 = [0:1:N/2-1,0, -N/2+1:-1]';
K12 = k1'.^2+k2.^2;
    
psi_hat              = fft2(w)./K12;
psi_hat(1,1,:)       = 0;
psi_hat(N/2+1,:,:)   = 0;
psi_hat(:,N/2+1,:)   = 0;
psi                  = real(ifft2(psi_hat));
u_hat =  1i*k2 .*psi_hat;
v_hat = -1i*k1'.*psi_hat;
u     = real(ifft2(u_hat));
v     = real(ifft2(v_hat));

end
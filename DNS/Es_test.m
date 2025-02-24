function [Es,kk] = Es_test(w,N)
%%
k1 = [0:N/2-1,0, -N/2+1:-1]';
k2 = [0:N/2-1,0, -N/2+1:-1]';
K12 = bsxfun(@(x,y) x.^2+y.^2,k1',k2);

w_hat            = fft2(w)/(N^2);
w_hat(1,1)       = 0;
w_hat(N/2+1,:)   = 0;
w_hat(:,N/2+1)   = 0;
psi_hat          = w_hat./K12;
psi_hat(1,1)     = 0;
psi_hat(N/2+1,:) = 0;
psi_hat(:,N/2+1) = 0;

Phik = 0.5* K12 .*abs(psi_hat).^2;
Phik = reshape(Phik,N^2,1);
%%
kk_max = round(sqrt(2)*(N/2-1));
kk = (1:kk_max)';
K12_index = reshape( sqrt(K12),[],1);
Es = zeros(length(kk),1);

for ind=1:length(kk)
    % index =  find( K12_index>=ind-0.5 & K12_index<ind+0.5 );
    Es(ind) = sum(Phik( K12_index>=ind-0.5 & K12_index<ind+0.5 ));
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
%
% compound trapezold integration, for (1/T)*\ind_{t0}^{t0+T} Ek(t) dt

function [Ek_mean] = meanEk_func(Ek)

N = length(Ek)-1;

% Ek_mean = sum(Ek(1:end-1) + Ek(2:end))/(2*N); % no faster

w = ones(N+1,1);w(1)=0.5;w(end)=0.5;
Ek_mean = sum(w.*Ek)/N;

end
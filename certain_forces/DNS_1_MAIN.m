%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%    Who am I?
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear;clc;
format long g;
%% physical space
N  = 2^7;
Lx = 2*pi;
Ly = 2*pi;
Re      = 400;
nu      = 1/Re;
[X,Y] = meshgrid((0:N-1)'*(Lx/N), (0:N-1)'*(Lx/N)); %(y,x)
%% spectral space
k1  = [0:1:N/2-1,0, -N/2+1:-1]';
k2  = [0:1:N/2-1,0, -N/2+1:-1]';
K12 = bsxfun(@(x,y) x.^2+y.^2 ,k1',k2);
%% initial conditions
%************************************************************************  
%******************  from decaying spectrum  **************************
%************************************************************************  
kf      = 3;
Ek_hat  = @(k1,k2) (9/11/kf)*( (sqrt(k1.^2+k2.^2)<=kf).*(k1.^2+k2.^2)./(kf^2)+ (sqrt(k1.^2+k2.^2)>kf).*((sqrt(k1.^2+k2.^2)/kf+eps).^(-5/3)) );
Ek_hat  = bsxfun(Ek_hat,[0:1:N/2, -N/2+1:-1],[0:1:N/2, -N/2+1:-1]')*N^2;
psi_hat = sqrt( 2*Ek_hat./(Lx*Ly*([0:1:N/2, -N/2+1:-1].^2+[0:1:N/2, -N/2+1:-1]'.^2 ) ) ).*exp( 1i*2*pi*rand(N,N) );

psi_hat(N/2+2:end,:)     = flipud(conj(psi_hat(2:N/2,:)));
psi_hat(N/2+2:end,2:end) = fliplr(psi_hat(N/2+2:end,2:end));
psi_hat(1,N/2+2:end)     = fliplr(conj(psi_hat(1,2:N/2)));
psi_hat(N/2+1,N/2+2:end) = fliplr(conj(psi_hat(N/2+1,2:N/2))); % the spectrum of real signals
psi_hat(1,1)             = 0;
psi_hat(N/2+1,N/2+1)     = 0; % (N/2,N/2)
psi_hat(:,N/2+1)         = 0; % (N/2,*)
psi_hat(N/2+1,:)         = 0; % (*,N/2)
w_hat                    = K12.*psi_hat;
w0                       = real(ifft2(w_hat));
psi0                     = real(ifft2(psi_hat));
clear('Ek_hat','w_hat','psi_hat','kf');

% external force
f12_hat  = f12hat_func(N);
%% Evolution process
parasF  = {Lx;Ly;N;nu};
T_force = 1;
dt      = 0.01;
t0      = 0;
t_max   = 20;
% [W_data,Ek,Dk,Ik,w,psi]                     = dns_2d_f_1(parasF,dt,t0,t_max,f12_hat,T_force,w0,psi0);
% [W_data,Psi_data,p_data,Ek,Dk,Ik,w,psi] = dns_2d_f_1_full(parasF,dt,t0,t_max,f12_hat,T_force,w0,psi0);
% [w,psi,Ek,Dk,Ik]                            = dns_2d_f_1_final(parasF,dt,t0,t_max,f12_hat,T_force,w0,psi0);
% [EE,Ek,Dk,Ik,w,psi]                         = dns_2d_f_1_proj(parasF,dt,t0,t_max,f12_hat,T_force,QQ,w0,psi0);
%% plotting 
close all;

COL1 = [linspace(0,1,2001)';1;ones(2001,1)];
COL2 = [linspace(0,1,2001)';1;linspace(1,0,2001)'];
COL3 = [ones(2001,1);1;linspace(1,0,2001)'];
MAP  = [COL1,COL2,COL3];
MAP(2002:2003,:) = [];
clear('COL1','COL2','COL3');
%% Saving AVI/MP4
t_start0 = 0;
Level    = 0.2;
fielName = 'W_test';
figure(4);
temp_size = get(0,'screensize');
set(gcf,'outerposition',[temp_size(1) temp_size(1) temp_size(4) temp_size(4)],...
        'color',[1 1 1],...
        'Toolbar','none');
videoMaker = VideoWriter(fielName,'MPEG-4');
open(videoMaker);
max_W = max(W_data(:));
min_W = min(W_data(:));
interval = 1;
if interval == 10
	formatSpec = '%.1f';
elseif interval==1
	formatSpec = '%.2f';
else
	formatSpec = '%.3f';
end
for ind=1:interval:size(W_data,3)
    [~,h] = contourf(X,Y,W_data(:,:,ind));
    h.LevelStep = Level;
    h.LineStyle = 'none';
    colormap(MAP);
    clim([min_W max_W]);
    colorbar;
    set(colorbar,'TickLength',0);
    set(gca,'LineWidth',2,...
            'FontName','Times New Roman',...
            'FontSize',22,...
            'TickLength',[0,0]);
    xlabel('$x$','Interpreter','latex');
    ylabel('$y$','Interpreter','latex');
    title(['$t=',num2str(t_start0+dt*(ind-1), formatSpec),'$'],'Interpreter','latex');
    axis equal;
    frame = getframe(gcf);
	writeVideo(videoMaker,frame);
end
close(videoMaker);
clear('videoMaker','t_start0','Level','fielName','temp_size','max_W','min_W','interval','h','frame','formatSpec');
close(4);
%% save as a movie, 'not visible'!!!
t_start0 = t1;
Level    = 0.1;
fielName = 'W_test';
fig1 = figure('visible','off');
temp_size = get(0,'screensize');
set(gcf,'outerposition',[temp_size(1) temp_size(1) temp_size(4) temp_size(4)],...
'color',[1 1 1],...
'Toolbar','none');
videoMaker = VideoWriter(fielName,'MPEG-4');
open(videoMaker);
max_W = max(W_data(:));
min_W = min(W_data(:));
interval = 100;
TOTAL = size(W_data,3);
p_old = 0;
for ind=1:interval:TOTAL
    [~,h] = contourf(X,Y,W_data(:,:,ind));
    h.LevelStep = Level;
    h.LineStyle = 'none';
    colormap(MAP);
    clim([min_W max_W]);
    colorbar;
    set(colorbar,'TickLength',0);
    set(gca,'LineWidth',2,...
    'FontName','Times New Roman',...
    'FontSize',22,...
    'TickLength',[0,0]);
    xlabel('$x$','Interpreter','latex');
    ylabel('$y$','Interpreter','latex');
    title(['$t=',num2str(t_start0+dt*(ind-1)),'$'],'Interpreter','latex');
    axis equal;
    frame = getframe(gcf);
    writeVideo(videoMaker,frame);
    p_new = round(ind/TOTAL*10);
    if p_new>p_old
        p_old = p_new;
        fprintf('------  %.1f%s has been done. ------\n',ind/TOTAL*100,'%');
    end
end
close(videoMaker);
close(fig1);
clear('videoMaker','t_start0','Level','fielName','temp_size','max_W','min_W','interval','h','frame','formatSpec','TOTAL','p_old','p_new');
clear('fig1')
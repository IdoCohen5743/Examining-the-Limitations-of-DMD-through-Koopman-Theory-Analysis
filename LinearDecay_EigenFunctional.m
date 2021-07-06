close all
clear all
clc
%% two nonorthogonal modes decay linearly 
%% signal
t = [0:0.1:50];
lambda1 = 1/30;
a1 = max((1-lambda1*t),0);
lambda2 = 1/10;
a2 = max((1-lambda2*t),0);
lengthy = 100;
V1 = zeros(lengthy,1);V1(20:70) = 1;
V2 = zeros(lengthy,1);V2(60:80) = 1;
VrefMode = [V2,V1];

Xall = V1*a1+V2*a2;
f = Xall(:,1);
%% data forming
M = length(t)-1;
X=Xall(:,1:M);Y=Xall(:,2:M+1);

%% SVD
r = 2;
[Ux,Sx,Vx] = svd(X);
figure();plot(diag(Sx))
%% dimentionality reduction
Uxr = Ux(:,[1:1:r]);
Vxr = Vx(:,[1:1:r]);
Sxr = Sx(1:1:r,1:1:r);
Cx = Uxr'*X;Cy = Uxr'*Y;
%% linear mapping
F = Cy*Cx'*inv(Cx*Cx');

[V,D,W] = eig(F);
D = diag(D);
phi = Uxr*V;
alpha = inv(V)*Uxr'*Xall(:,1);
aphi = phi*diag(alpha);
figure();plot(V2);hold on;plot(real(aphi(:,1)));hold on;plot(imag(aphi(:,1)));
figure();plot(V1);hold on;plot(real(aphi(:,2)));hold on;plot(imag(aphi(:,2)));
%% ref 
h1 = figure();movegui('northwest');elc = get(gca,'colororder');
plot([-lengthy/2:1:lengthy/2-1],VrefMode(:,1),'Color',elc(2,:),'LineStyle','--','LineWidth',8);
set(gca,'TickLabelInterpreter','latex');
set(gca,'FontSize',52);grid on;
xlim([-lengthy/2,lengthy/2])
set(gca,'XTick',[-50:25:50]);
set(gca,'YLim',[-0.1,1.2]);
set(gca,'YTick',[-0.8:0.4:0.8]);
ylabel('$v_1(x)$','Interpreter','latex','Position', [-lengthy*6/10 0.95],'Rotation',0);%,'FontSize',95);
xlabel('$x$','Interpreter','latex','Position', [lengthy*5.5/10, 0]);%,'FontSize',95);
pos = [0.15    0.1100    0.7750    0.8150];
set(gca,'Position',pos)
write_pdf_New_Image(['RefMode1'],h1,1000,750)

h1 = figure();movegui('northwest');elc = get(gca,'colororder');
plot([-lengthy/2:1:lengthy/2-1],VrefMode(:,2),'Color',elc(2,:),'LineStyle','--','LineWidth',8);
set(gca,'TickLabelInterpreter','latex');
set(gca,'FontSize',52);grid on;
xlim([-lengthy/2,lengthy/2])
set(gca,'XTick',[-50:25:50]);
set(gca,'YLim',[-0.1,1.2]);
set(gca,'YTick',[-0.8:0.4:0.8]);
ylabel('$v_2(x)$','Interpreter','latex','Position', [-lengthy*6/10 0.95],'Rotation',0);%,'FontSize',95);
xlabel('$x$','Interpreter','latex','Position', [lengthy*5.5/10, 0]);%,'FontSize',95);
pos = [0.15    0.1100    0.7750    0.8150];
set(gca,'Position',pos)
write_pdf_New_Image(['RefMode2'],h1,1000,750)



%% display DMD mode
h1 = figure();movegui('northwest');
% subplot(1,5,1);plot([-lengthy/2:1:lengthy/2-1],sum(aphi,2),'LineWidth',8);hold on;plot([-lengthy/2:1:lengthy/2-1],f,'--','LineWidth',8);hold off;
% title(['$\,$'],'Interpreter', 'latex')
% set(gca,'TickLabelInterpreter','latex')
% set(gca,'FontSize',35);grid on;
% xlim([-lengthy/2,lengthy/2])
% set(gca,'XTick',[-50:25:50]);
% ylabel('$u(x,t=0)$','Interpreter','latex','Position', [-lengthy*6/10 2.2],'Rotation',0);%,'FontSize',95);
% xlabel('$x$','Interpreter','latex','Position', [lengthy*6/10, -0.37]);%,'FontSize',95);
% set(gca,'YLim',[-0.5,2.2]);
% set(gca,'YTick',[-0.4:0.4:2.3]);
for iii=2:2:4
    subplot(1,4,iii-1);
    plot([-lengthy/2:1:lengthy/2-1],real(aphi(:,iii/2)),'LineWidth',8);hold on;plot([-lengthy/2:1:lengthy/2-1],VrefMode(:,iii/2),'--','LineWidth',8);hold on
%     title(['$Re\{\phi_',num2str(iii/2),'\}$'],'Interpreter', 'latex')
    ylabel('$Re$','Interpreter','latex','Position', [-lengthy*7/10 0.9],'Rotation',0);%,'FontSize',95);
    set(gca,'TickLabelInterpreter','latex')
    set(gca,'FontSize',35);grid on;
    xlim([-lengthy/2,lengthy/2])
    set(gca,'XTick',[-50:25:50]);
    set(gca,'YLim',[-1.2,1.2]);
    set(gca,'YTick',[-0.8:0.4:0.8]);
    
    subplot(1,4,iii);plot([-lengthy/2:1:lengthy/2-1],imag(aphi(:,iii/2)),'LineWidth',8);
%     title(['$Im\{\phi_',num2str(iii/2),'\}$'],'Interpreter', 'latex')
    ylabel('$Im$','Interpreter','latex','Position', [-lengthy*7/10 0.9],'Rotation',0);%,'FontSize',95);
    set(gca,'TickLabelInterpreter','latex')
    set(gca,'FontSize',35);grid on;
    xlim([-lengthy/2,lengthy/2])
    set(gca,'XTick',[-50:25:50]);
    set(gca,'YLim',[-1.2,1.2]);
    set(gca,'YTick',[-0.8:0.4:0.8]);
    
end
write_pdf_New_Image(['modesDMDLinear'],h1,2670,430)



%% flow reconstruction
Xkova = phi(:,1)*alpha(1)*D(1).^[0:1:M]+phi(:,2)*alpha(2)*D(2).^[0:1:M];
% for iii=1:1:M
%     h=figure(26);
%     plot(Xall(:,iii)); hold on; plot(real(Xkova(:,iii)));hold off;set(gca,'YLim',[-0.2000 2.2]);
%     drawnow;
% end

%% display reconstruction
[X,Y] = meshgrid(t,[-lengthy/2:1:lengthy/2-1]);
h1 = figure();movegui('northwest');surf(X,Y,Xall,'LineStyle','none');
set(gca,'CameraPosition',[365.4171 -314.5867    9.6603]);
zlabel('$u(x,t)$','Interpreter','latex','Position', [-3.0241  -54.9597    2.15],'Rotation',0);
ylabel('$x$','Interpreter','latex','Position', [57.0465   25.0434   -0.0793]);
xlabel('$t$','Interpreter','latex','Position', [24.4680  -57.3709   -0.0736]);

set(gca,'TickLabelInterpreter','latex')
set(gca,'XLim',[0,t(end)]);
set(gca,'FontSize',52);
write_pdf_New_Image(['Xall'],h1,1000,750)

h1 = figure();movegui('northwest');surf(X,Y,Xkova,'LineStyle','none');
set(gca,'CameraPosition',[365.4171 -314.5867    9.6603]);
zlabel('$\hat{u}(x,t)$','Interpreter','latex','Position', [-3.0241  -54.9597    2.15],'Rotation',0);
ylabel('$x$','Interpreter','latex','Position', [57.0465   25.0434   -0.0793]);
xlabel('$t$','Interpreter','latex','Position', [24.4680  -57.3709   -0.0736]);

set(gca,'TickLabelInterpreter','latex')
set(gca,'XLim',[0,t(end)]);
set(gca,'FontSize',52);
write_pdf_New_Image(['Xkova'],h1,1000,750)

h1 = figure();movegui('northwest');surf(X,Y,Xkova-Xall,'LineStyle','none');
set(gca,'CameraPosition',[446.8843 -344.9111    1.2052]);
ylabel('$x$','Interpreter','latex','Position', [57.0465   25.0434   -0.0793]);
xlabel('$t$','Interpreter','latex','Position', [24.4680  -57.3709   -0.0736]);

set(gca,'TickLabelInterpreter','latex')
set(gca,'XLim',[0,t(end)]);
set(gca,'FontSize',52);
write_pdf_New_Image(['Xdiff'],h1,1000,750)
%% 
%% decay dictionary
%% sparse representation
% init
A = [];
for iii=1:1:40
    lambda = 1/iii;
    aTemp = max((1-lambda*t),0);
    A = [A;aTemp];
end

param.lambda=0.15; % not more than 20 non-zeros coefficients
param.numThreads=-1; % number of processors/cores to use; the default choice is -1
                     % and uses all the cores of the machine
param.mode=2;        % penalized formulation

V=full(mexLasso(Xall',A',param)');
[~,col]=size(V);
%% remove the zero modes
iii=1;
while(true)
    [~,col]=size(V);
    if(iii>col)
        break;
    end
    
    if(norm(V(:,iii))<1e-10)
        V(:,iii)=[];
        A(iii,:)=[];
    end
    iii = iii+1;
end
Ainit = A;
Vinit = V;
%%
Err = [];
while(true) 
    Err = [Err,norm(V*A-Xall,'f')];
    rating = vecnorm(V);
    [values,ind] = sort(rating);
%     V(:,ind(1))=[];
    A(ind(1),:)=[];
    [row,~] = size(A);
    row
    if(row<1)
        break;
    end
    V = Xall*A'*inv(A*A');
end

% init
A = Ainit;
V = Vinit;
while(true) 
    if(norm(V*A-Xall,'f')<=min(Err))
        break;
    end
    rating = vecnorm(V);
    [values,ind] = sort(rating);
%     V(:,ind(1))=[];
    A(ind(1),:)=[];
    [row,~] = size(A);
    row
    if(row<1)
        break;
    end
    V = Xall*A'*inv(A*A');
end
figure();plot(Err)
figure();plot(A(1,:));hold on;plot(a2)
figure();plot(A(2,:));hold on;plot(a1)

figure();plot(V(:,1));hold on; plot(V2);
figure();plot(V(:,2));hold on; plot(V1);
%% display
h1 = figure();movegui('northwest');
% subplot(1,3,1);plot([-lengthy/2:1:lengthy/2-1],sum(aphi,2),'LineWidth',8);hold on;plot([-lengthy/2:1:lengthy/2-1],f,'--','LineWidth',8);hold off;
% set(gca,'TickLabelInterpreter','latex')
% set(gca,'FontSize',35);grid on;
% xlim([-lengthy/2,lengthy/2])
% set(gca,'XTick',[-50:25:50]);
% ylabel('$u(x,t=0)$','Interpreter','latex','Position', [-lengthy*5.5/10 2.2],'Rotation',0);%,'FontSize',95);
% xlabel('$x$','Interpreter','latex','Position', [lengthy*6/10, 0]);%,'FontSize',95);
% set(gca,'YLim',[-0.1,2.2]);
% set(gca,'YTick',[-0.4:0.4:2.2]);
for iii=2:1:3
    subplot(1,2,iii-1);
    plot([-lengthy/2:1:lengthy/2-1],V(:,iii-1),'LineWidth',8);hold on;plot([-lengthy/2:1:lengthy/2-1],VrefMode(:,iii-1),'--','LineWidth',8);hold on
%     title(['$\phi_',num2str(iii-1),'$'],'Interpreter', 'latex')
ylabel(['$v_',num2str(iii-1),'(x)$'],'Interpreter','latex','Position', [-lengthy*5.5/10 1.1],'Rotation',0);%,'FontSize',95);
    set(gca,'TickLabelInterpreter','latex')
    set(gca,'FontSize',35);grid on;
    xlim([-lengthy/2,lengthy/2])
    set(gca,'XTick',[-50:25:50]);
    set(gca,'YLim',[-0.1,1.1]);
    set(gca,'YTick',[-0.4:0.2:2.2]);
end
title(['$\, $'],'Interpreter', 'latex')
write_pdf_New_Image(['sparseLinear'],h1,1068,430)

%% flow reconstruction
Xkova = V*A;
% for iii=1:1:M
%     h=figure(26);
%     plot(Xall(:,iii)); hold on; plot(real(Xkova(:,iii)));hold off;set(gca,'YLim',[-0.2000 2.2]);
%     drawnow;
% end
%% Display reconstruction
[X,Y] = meshgrid(t,[-lengthy/2:1:lengthy/2-1]);
% h1 = figure();movegui('northwest');surf(X,Y,Xall,'LineStyle','none');
% set(gca,'CameraPosition',[365.4171 -314.5867    9.6603]);
% zlabel('$y$','Interpreter','latex','Position', [-3.0241  -54.9597    1.1181],'Rotation',0);
% ylabel('$x$','Interpreter','latex','Position', [57.0465   25.0434   -0.0793]);
% xlabel('$t$','Interpreter','latex','Position', [24.4680  -57.3709   -0.0736]);
% 
% set(gca,'TickLabelInterpreter','latex')
% set(gca,'XLim',[0,t(end)]);
% set(gca,'FontSize',52);
% write_pdf_New_Image(['XallSparse'],h1,1000,750)

h1 = figure();movegui('northwest');surf(X,Y,Xkova,'LineStyle','none');
set(gca,'CameraPosition',[365.4171 -314.5867    9.6603]);
zlabel('$\hat{u}(x,t)$','Interpreter','latex','Position', [-3.0241  -54.9597    2.15],'Rotation',0);
ylabel('$x$','Interpreter','latex','Position', [57.0465   25.0434   -0.0793]);
xlabel('$t$','Interpreter','latex','Position', [24.4680  -57.3709   -0.0736]);

set(gca,'TickLabelInterpreter','latex')
set(gca,'XLim',[0,t(end)]);
set(gca,'FontSize',52);
write_pdf_New_Image(['XkovaSparse'],h1,1000,750)

h1 = figure();movegui('northwest');surf(X,Y,Xkova-Xall,'LineStyle','none');
set(gca,'CameraPosition',[387.0996 -421.9767    0.0000]);
zlabel('$u(x,t)$','Interpreter','latex','Position', [-3.0241  -54.9597    1.1181],'Rotation',0);
ylabel('$x$','Interpreter','latex','Position', [57.0465   25.0434   -0.0793]);
xlabel('$t$','Interpreter','latex','Position', [24.4680  -57.3709   -0.0736]);

set(gca,'TickLabelInterpreter','latex')
set(gca,'XLim',[0,t(end)]);
set(gca,'FontSize',52);
write_pdf_New_Image(['XdiffSparse'],h1,1000,750)


%% eigenfunctionals
trecon = -(inv(V'*V)*V'*Xall-1);

J1 = exp(trecon(1,:));J2 = exp(trecon(2,:));
h1 = figure();movegui('northwest');
plot(t,trecon(1,:),'LineWidth',8);
set(gca,'YLim',[0,1]);
pos = get(gca,'Position');
pos(1) = 0.15;
set(gca,'Position',pos);
set(gca,'FontSize',80);grid on;
xlim([0,t(end)])
set(gca,'TickLabelInterpreter','latex')
ylabel('$\ln(\phi_1)$','Interpreter','latex','Position', [-5 0.825],'FontSize',95,'Rotation',0);
xlabel('$t$','Interpreter','latex','Position', [t(end)+1 0.075],'FontSize',95);
write_pdf_New_Image(['J1'],h1,2000,1000)


h1=figure();movegui('northwest');
plot(t,trecon(2,:),'LineWidth',8);
set(gca,'YLim',[0,1]);
pos = get(gca,'Position');
pos(1) = 0.15;
set(gca,'Position',pos);
set(gca,'FontSize',80);grid on;
xlim([0,t(end)])
set(gca,'TickLabelInterpreter','latex')
ylabel('$\ln(\phi_2)$','Interpreter','latex','Position', [-5 0.825],'FontSize',95,'Rotation',0);
xlabel('$t$','Interpreter','latex','Position', [t(end)+1 0.075],'FontSize',95);
write_pdf_New_Image(['J2'],h1,2000,1000)
return

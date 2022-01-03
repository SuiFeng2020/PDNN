% Computing the error for different parameter 
load('Gmsh_262_cos_Snapshots_time_parameter05.mat')
load('Gmsh_262_cos_PODbasis_twostep_3time_4parameter05.mat')
load('Gmsh_262_cos_tptrain05.mat')
load dof_mesh_26_2.mat;
load Gmsh_262_dgtdmatrix.mat
%% for a test parameter
% %% DGTD solutions
% global parameter;
% %epsilon_t and epsilon_mu at test prameter epsilon = 2.215
% epsilon0 = 8.854*(1e-12); 
% mu0 = 4*pi*(10^(-7));
% parameter.c0 = sqrt(1/(epsilon0*mu0));
% parameter.c = 1;
% parameter.freq = 3e8; % Hz  
% parameter.dt = train.dt;
% parameter.tmax = 50/parameter.freq; % the max time;
% parameter.pOrder = 2; % interpolation order of DGTD
% parameter.prob = 3; % mesh type
% %
% % testparater
% % 1th ;
% test.parameter = 2.215;
% parameter.re.epstest = test.parameter;
% %
% Nttr = length(train.time); 
% Nutr = length(train.parameter);
% % call the main subroutine
% DGTDtime = DGTDmultilayer;
% %% get projection error 
% kk = 0; % get the snapshots with testesp = 1:0.05:5;
% epsilon_t = [1e-1 1e-2 1e-3 1e-4];
% epsilon_mu = [1e-1 1e-2 1e-3 1e-4];
% for iet = 1:length(epsilon_t)
%     % the 1th svd
%     for ii = 1:Nutr
%         % computing the POD basis for the ii parameter
%         tPOD(ii) = getPODbasis(1 - epsilon_t(iet),...
%                                Snapshots(ii)); 
%         ii
%     end
%     cattPOD.Basis.Hxe = [];
%     cattPOD.Basis.Hye = [];
%     cattPOD.Basis.Eze = [];
%     for ii = 1:Nutr
%         cattPOD.Basis.Hxe = cat(2,cattPOD.Basis.Hxe,tPOD(ii).Basis.Hx);
%         cattPOD.Basis.Hye = cat(2,cattPOD.Basis.Hye,tPOD(ii).Basis.Hy);
%         cattPOD.Basis.Eze = cat(2,cattPOD.Basis.Eze,tPOD(ii).Basis.Ez);
%     end
%     %
%     for iemu = 1:length(epsilon_mu)
%         % the 2th svd
%         timeparameterPOD = getPODbasis(1 - epsilon_mu(iemu),cattPOD.Basis);
%         protimeErrorL2 = zeros(1,Nttr);
%         for jj  = 1:Nttr
%             DGTDTime = [DGTDtime.Hxe(:,jj),DGTDtime.Hye(:,jj),DGTDtime.Eze(:,jj)];
%             proMORCSITime = [timeparameterPOD.Basis.Hx*(timeparameterPOD.Basis.Hx'*DGTDtime.Hxe(:,jj)),...
%                              timeparameterPOD.Basis.Hy*(timeparameterPOD.Basis.Hy'*DGTDtime.Hye(:,jj)),...
%                              timeparameterPOD.Basis.Ez*(timeparameterPOD.Basis.Ez'*DGTDtime.Eze(:,jj))];
%             [proerrE, proerrH] = getErr(proMORCSITime,DGTDTime);
%             error.protimeErrorL2(1,jj) = sqrt(proerrE^2 + proerrH^2);
%         end
%         apProerrorsvd(iet,iemu) = sum(error.protimeErrorL2)/Nttr;
%     end    
% end
% % 
% save proerror.mat epsilon_t epsilon_mu apProerrorsvd
% apProerrorsvd = Proerrorsvd;
% figure(1)
% h = gca;
% set(h,'FontSize',20);
% [mt,mmu] = meshgrid(epsilon_t,epsilon_mu);
% [x,y] = meshgrid(linspace(0.0001,0.1,500),linspace(0.0001,0.1,500));
% [m,n] = size(apProerrorsvd);
% reshmt = reshape(mt',m*n,1);
% reshmmu = reshape(mmu',m*n,1);
% reshPro = reshape(apProerrorsvd,m*n,1);
% fPro = TriScatteredInterp(reshmt,reshmmu,reshPro);
% surf(x,y,fPro(x,y))
% % plot(0.001,0.0001,'r*','MarkerSize',30)
% % [mt,mmu] = meshgrid(epsilon_t,epsilon_mu);
% % surf(mt',mmu',Proerrorsvd)
% set(gca,'xscale','log')
% set(gca,'yscale','log')
% shading interp; 
% axis tight
% colorbar; 
% view(2);
% h = colorbar;
% set(h,'fontsize',20);
% t1 = caxis;
% t1 = linspace(t1(1),t1(2),4);
% my_handle = colorbar('ytick',t1,'fontsize',20);
% xlabel('$$\rho_t$$','Interpreter','Latex')
% ylabel('$$\rho_\theta$$','Interpreter','Latex')
% set(gcf,'PaperPositionmode','auto')
% print('-depsc','cylindergmsh3cossplineproerrortwosvdap2215')
%%  for all train parameter
Nttr = length(train.time); 
Nutr = length(train.parameter);
%% get projection error 
kk = 0; % get the snapshots with testesp = 1:0.05:5;
epsilon_t = [1e-1 1e-2 1e-3 1e-4];
epsilon_mu = [1e-1 1e-2 1e-3 1e-4];
allpProerrorsvd = zeros(length(epsilon_t),length(epsilon_mu));
allpProerrorsvdE = zeros(length(epsilon_t),length(epsilon_mu));
allpProerrorsvdH = zeros(length(epsilon_t),length(epsilon_mu));
Ndof = size(ADGTD.Me,1);
zeronDGTDTime = zeros(Ndof,3);
for iet = 1:length(epsilon_t)
    %% the 1th svd
    for ii = 1:Nutr
        % computing the POD basis for the ii parameter
        tPOD(ii) = getPODbasis(1 - epsilon_t(iet),...
                               Snapshots(ii)); 
        ii
    end
    cattPOD.Basis.Hxe = [];
    cattPOD.Basis.Hye = [];
    cattPOD.Basis.Eze = [];
    for ii = 1:Nutr
        cattPOD.Basis.Hxe = cat(2,cattPOD.Basis.Hxe,tPOD(ii).Basis.Hx);
        cattPOD.Basis.Hye = cat(2,cattPOD.Basis.Hye,tPOD(ii).Basis.Hy);
        cattPOD.Basis.Eze = cat(2,cattPOD.Basis.Eze,tPOD(ii).Basis.Ez);
    end
    %%
    for iemu = 1:length(epsilon_mu)
        % the 2th svd
        timeparameterPOD = getPODbasis(1 - epsilon_mu(iemu),cattPOD.Basis);
%         error.protimeErrorL2 = zeros(Nutr,Nttr);
        error.reprotimeErrorE = zeros(Nutr,Nttr);
        error.reprotimeErrorH = zeros(Nutr,Nttr);
        for is = 1:Nutr
            DGTDtime = Snapshots(is); % for is-th parameter
            for jj  = 1:Nttr
                 DGTDTime = [DGTDtime.Hxe(:,jj),DGTDtime.Hye(:,jj),DGTDtime.Eze(:,jj)];
                 proMORCSITime = [timeparameterPOD.Basis.Hx*(timeparameterPOD.Basis.Hx'*DGTDtime.Hxe(:,jj)),...
                                  timeparameterPOD.Basis.Hy*(timeparameterPOD.Basis.Hy'*DGTDtime.Hye(:,jj)),...
                                  timeparameterPOD.Basis.Ez*(timeparameterPOD.Basis.Ez'*DGTDtime.Eze(:,jj))];
                 [proerrE, proerrH] = getErr(proMORCSITime,DGTDTime);
                 [reproerrE, reproerrH] = getErr(zeronDGTDTime,DGTDTime);
%                  error.protimeErrorL2(is,jj) = sqrt(proerrE^2 + proerrH^2);
                 error.reprotimeErrorE(is,jj) = proerrE/reproerrE;
                 error.reprotimeErrorH(is,jj) = proerrH/reproerrH;
            end
        end
%         allpProerrorsvd(iet,iemu) = sum(sum(error.protimeErrorL2,2))/(Nttr*Nutr);  
        allpProerrorsvdE(iet,iemu) = sum(sum(error.reprotimeErrorE,2))/(Nttr*Nutr);
        allpProerrorsvdH(iet,iemu) = sum(sum(error.reprotimeErrorH,2))/(Nttr*Nutr);
    end    
end
% 
save reallpproerrorEH.mat epsilon_t epsilon_mu allpProerrorsvdE allpProerrorsvdH
%% show the projection error
figure(2)
h = gca;
set(h,'FontSize',20);
[mt,mmu] = meshgrid(epsilon_t,epsilon_mu);
[x,y] = meshgrid(linspace(0.0001,0.1,500),linspace(0.0001,0.1,500));
[m,n] = size(allpProerrorsvdE);
reshmt = reshape(mt',m*n,1);
reshmmu = reshape(mmu',m*n,1);
reshPro = reshape(allpProerrorsvdE,m*n,1);
fPro = TriScatteredInterp(reshmt,reshmmu,reshPro);
surf(x,y,fPro(x,y))
% plot(0.001,0.0001,'r*','MarkerSize',30)
% [mt,mmu] = meshgrid(epsilon_t,epsilon_mu);
% surf(mt',mmu',Proerrorsvd)
set(gca,'xscale','log')
set(gca,'yscale','log')
shading interp; 
axis tight
colorbar; 
view(2);
h = colorbar;
set(h,'fontsize',20);
t1 = caxis;
t1 = linspace(t1(1),t1(2),4);
my_handle = colorbar('ytick',t1,'fontsize',20);
xlabel('$$\varrho_t$$','Interpreter','Latex')
ylabel('$$\varrho_\theta$$','Interpreter','Latex')
set(gcf,'PaperPositionmode','auto')
print('-depsc','cylindergmsh3cossplinereproerrorEtwosvdallp')
%
figure(3)
h = gca;
set(h,'FontSize',20);
[mt,mmu] = meshgrid(epsilon_t,epsilon_mu);
[x,y] = meshgrid(linspace(0.0001,0.1,500),linspace(0.0001,0.1,500));
[m,n] = size(allpProerrorsvdH);
reshmt = reshape(mt',m*n,1);
reshmmu = reshape(mmu',m*n,1);
reshPro = reshape(allpProerrorsvdH,m*n,1);
fPro = TriScatteredInterp(reshmt,reshmmu,reshPro);
surf(x,y,fPro(x,y))
% plot(0.001,0.0001,'r*','MarkerSize',30)
% [mt,mmu] = meshgrid(epsilon_t,epsilon_mu);
% surf(mt',mmu',Proerrorsvd)
set(gca,'xscale','log')
set(gca,'yscale','log')
shading interp; 
axis tight
colorbar; 
view(2);
h = colorbar;
set(h,'fontsize',20);
t1 = caxis;
t1 = linspace(t1(1),t1(2),4);
my_handle = colorbar('ytick',t1,'fontsize',20);
xlabel('$$\varrho_t$$','Interpreter','Latex')
ylabel('$$\varrho_\theta$$','Interpreter','Latex')
set(gcf,'PaperPositionmode','auto')
print('-depsc','cylindergmsh3cossplinereproerrorHtwosvdallp')
   


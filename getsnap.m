function Snapshots = getsnap(eps,parameter)
%
global refTri dofmat ADGTD;
%
nTranges = size(dofmat.DOF(:,:,1),1);
% Initialization incident wave
eInc = zeros(1,nTranges*refTri.nde);
hInc = zeros(2,nTranges*refTri.nde);
omega = 2*pi*parameter.freq/parameter.c0; % normalied
k = sqrt(1*1)*omega/(parameter.c); % the normalized speed of c = 1;
Hyinc = @(x,y,t)(-cos(omega*t - k*x));
Ezinc = @(x,y,t)(cos(omega*t - k*x));
% Hyinc = @(x,y,t)(-exp((1i)*omega*t - (1i)*k*x));
% Ezinc = @(x,y,t)(exp((1i)*omega*t - (1i)*k*x));
% Initialization incident wave
xDof = reshape(dofmat.DOF(:,:,1)',nTranges*refTri.nde,1);
yDof = reshape(dofmat.DOF(:,:,2)',nTranges*refTri.nde,1);
% Initial conditions
Hxe = zeros(nTranges*refTri.nde,1);
Hye = zeros(nTranges*refTri.nde,1);
Eze = zeros(nTranges*refTri.nde,1);
%% get the MOR matrix except ADGTD.Me
MetalDof = reshape(ADGTD.metaldof,refTri.nde*size(ADGTD.metaldof,2),1);
ADGTDMeeps = ADGTD.Me;
ADGTDMeeps(MetalDof,MetalDof) = ADGTDMeeps(MetalDof,MetalDof)*eps;
%%
tic
tt = 0;
time = 0;
deltat = parameter.dt;
ell = 1;
Ntmax = parameter.tmax*parameter.freq - 1;
% dgtdtime = 0:deltat:parameter.tmax*parameter.c0;
% critime = dgtdtime(size(dgtdtime,2) - size(find(dgtdtime>=Ntmax),2) + 1);
%
while time <= parameter.tmax*parameter.c0
    %
    eInc(1,dofmat.vmapB) = Ezinc(xDof(dofmat.vmapB),yDof(dofmat.vmapB),...
        deltat*(tt));
    hInc(2,dofmat.vmapB) = Hyinc(xDof(dofmat.vmapB),yDof(dofmat.vmapB),...
        deltat*(tt+0.5));
    Eze = Eze + ((ADGTDMeeps)\...
        (0.5*ADGTD.Ky*Hxe - 0.5*ADGTD.SyI*Hxe - 0.5*(ADGTD.SyB.y*Eze ...
        + ADGTD.SyBB*hInc(1,:)' -  ADGTD.SyB.y*eInc') ...
        - 0.5*ADGTD.Kx*Hye + 0.5*ADGTD.SxI*Hye + 0.5*(-ADGTD.SxB.x*Eze ...
        + ADGTD.SxBB*hInc(2,:)' + ADGTD.SxB.x*eInc')))*deltat;
    %    eInc(1,dofmat.vmapB) = Ezinc(xDof(dofmat.vmapB),yDof(dofmat.vmapB),...
    %                                                    deltat*(tt + 1));
    Hxe = Hxe + ((ADGTD.Mh)\...
        (0.5*ADGTD.Ky*Eze - 0.5*ADGTD.SyI*Eze - 0.5*((ADGTD.SyB.y*Hxe ...
        - ADGTD.SyB.x*Hye) + ADGTD.SyBB*eInc' - (ADGTD.SyB.y*hInc(1,:)' ...
        - ADGTD.SyB.x*hInc(2,:)'))))*deltat;
    Hye = Hye + ((ADGTD.Mh)\...
        (-0.5*ADGTD.Kx*Eze + 0.5*ADGTD.SxI*Eze + 0.5*((ADGTD.SxB.y*Hxe ...
        - ADGTD.SxB.x*Hye) + ADGTD.SxBB*eInc' - (ADGTD.SxB.y*hInc(1,:)' ...
        - ADGTD.SxB.x*hInc(2,:)'))))*deltat;
    %     dgsolution = [Hxe, Hye, Eze];
    %   Increment time
    time  = time + deltat;
    tt = tt + 1;
    %
    %     if abs(time - critime) < 1e-10;
    %         Snapshots.Hxe(:,ell) = Hxe;
    %         Snapshots.Hye(:,ell) = Hye;
    %         Snapshots.Eze(:,ell) = Eze;
    %         Snapshots.traintime(1,ell) = time;
    %         critime = critime + parameter.snapell*deltat;
    %         ell = ell + 1;
    %    end
    if time >= Ntmax
        Snapshots.Hxe(:,ell) = Hxe;
        Snapshots.Hye(:,ell) = Hye;
        Snapshots.Eze(:,ell) = Eze;
        Snapshots.traintime(1,ell) = time;
        ell = ell + 1;
    end
    
end
toc

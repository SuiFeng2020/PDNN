function [Snapshots,tPOD,timeparameterPOD,train] = DGTDPODSolutionCylinder(parameter)
% Input : pOrder  -- interpolation order of DGTD
%         visual  -- visualization parameter: 
%                    0 -> plot nothing; 
%                    1 -> plot the field value on the x-axis
%                         and the y-axis for Ex, Ez and Hy; 
%                    2 -> plot the intensity distribution
%                         of Ex, Ez and Hy; 
%                    3 -> plot both 1 and 2
%         prob    -- problem type (mesh type):
%                    0 -> coarse mesh
%                    1 -> refined mesh (once)
%         tmax    -- the final time
% Output: dataset:
%         ell     -- the number of snapshots for any parameter
%         Au      -- the all snapshots
%         POD     -- the POD bases
% Author: K. Li
% Date  : 2019-7-11
%--------------------------------------------------
%% Global variable
global totMsh;
global totEdg; 
global refTri;
global dofmat;
global ADGTD;
%% Set up
t1 = clock;
% meshing
tic
disp('Mesh creating ...')
% mesh loading
totMsh = getMesh(parameter); % gmsh
% totMsh = gTriMesh(parameter); % matlab 
% find out the faces
totEdg = getFace;
% get the information of matrix (mass, stiff,..) in reference element
refTri = info_Pk(parameter); 
% get the dof
dofmat = getdof;  
disp('CPU time for mesh creating:')
tMesh = toc
%%  The matrices of assembly
tic
disp('Matrix assembling ...')                 
ADGTD = DGTDPODAssembly;
disp('CPU time for matrix assembly:') 
tAssembly = toc
%% Compute time step size
CFL = [0.3 0.2 0.1 0.08 0.06];%CFL-value
epsilon0 = 8.854*(1e-12);
mu0 = 4*pi*(10^(-7));
parameter.c0 = sqrt(1/(epsilon0*mu0));
parameter.c = 1;
% betaik = (parameter.pOrder + 1)*(parameter.pOrder + 2)/2;
al_be = (8+sqrt(40))/3;
lenpara = length(parameter.re.epstest);
train.parameter = parameter.re.epstest;
dt = zeros(1,lenpara);
for ii = 1:lenpara
    eps = train.parameter(ii);
    dtscale = getdt(eps);
    dt(1,ii) = 0.9*CFL(parameter.pOrder)*(4*min(dtscale));
end
train.dt = min(dt);
parameter.dt = train.dt; % the delta t of DGTD method for all parameter
disp('CPU time for construction:')
t2 = clock;
etime(t2,t1)
%% getting snapshots for the ii parameter
for ii = 1:lenpara
    eps = train.parameter(ii);
    Snapshots(ii) = getsnap(eps,parameter);
%     tPOD(ii) = getPODbasis(parameter.info1,...
%            Snapshots(ii)); % computing the POD basis for the ii parameter
%     ii
end
train.time = Snapshots(1).traintime;
% We have obtained the all snapshots
tic
kk = 0; % get the snapshots with testesp = 1:0.05:5;
for ii = 1:lenpara
    tPOD(ii) = getPODbasis(1 - 1e-3,...
           Snapshots(ii+kk)); % computing the POD basis for the ii parameter
%     kk = kk +1
    ii
end
%% getting the POD basis based on the twostep SVD for all papameter
cattPOD.Basis.Hxe = [];
cattPOD.Basis.Hye = [];
cattPOD.Basis.Eze = [];
for ii = 1:lenpara
    cattPOD.Basis.Hxe = cat(2,cattPOD.Basis.Hxe,tPOD(ii).Basis.Hx);
    cattPOD.Basis.Hye = cat(2,cattPOD.Basis.Hye,tPOD(ii).Basis.Hy);
    cattPOD.Basis.Eze = cat(2,cattPOD.Basis.Eze,tPOD(ii).Basis.Ez);
end
timeparameterPOD = getPODbasis(1 - 1e-4,cattPOD.Basis);
toc
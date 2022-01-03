function dgsolution = DGTDSolutionCylinder
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
global parameter;
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
disp('CPU time for construction:')
t2 = clock;
etime(t2,t1)
%%
eps = parameter.re.epstest;
dgsolution = getsnap(eps,parameter); % get snapshots


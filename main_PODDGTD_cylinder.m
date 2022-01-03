%% main_PODDGTD_cylider
% This program simulates Scattering of plane wave by 
% a dielectric cylinder
% The radius of the disk = 0.6 m;
% The radius of the AAB (artificial absorbing boundary): 
% [-1.6 1.6]\times[-1.6 1.6]
% Author: K.Li        
% Date  : 2019-7-11 
% Modified : 2019-7-11
%%    
% reset Matlab workspace     
 clear; clc; close all;  
% set parameters    
freq = 3e8; % Hz       
tmax = 50/freq; % the max time;
pOrder = 2; % interpolation order of DGTD
prob = 3; % mesh type
visual = 2; % show the mesh plot,  in Visdgtdsolution.m
fixed_point = 1; 
re.epstest = 1:0.05:5;
%re.epstest = [1.215, 2.215, 3.215, 4.215];
info1 = 1 - 1e-3; % uses-specifed tolerance for time SVD
info2 = 1 - 1e-3; % uses-specifed tolerance for parameter SVD
% call the main subroutine
parameter = struct('prob',prob,'freq',freq,...
                   'pOrder',pOrder,'re',re,...
                   'tmax',tmax,'visual',visual,...
                   'fixed_point',fixed_point,...
                   'info1',info1,'info2',info2);
[Snapshots,tPOD, timeparameterPOD,train] = DGTDPODSolutionCylinder(parameter);
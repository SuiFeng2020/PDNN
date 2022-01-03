function Pk_ref = info_Pk(parameter)

% Function giving the parameters of the Pk family on the reference
% triangle.
% USAGE:
% Pk_ref = info_Pk(p)
%
% -------------------------------------------------------------------------
% Reference element.
% 3 (0, 1)
% |`
% |  `
% |    `
% |      `
% |        `
% 1----------2
% (0, 0)     (1, 0)
% Edges orientation. edge 1: 1->2, edge 2: 2->3, edge 3: 3->1.
% -------------------------------------------------------------------------
% INPUT.
% p: polynomial order.
%
% OUTPUT.
% Pk_ref: structure for the reference element and containing
% * nde: Number of degrees of freedom (dof) in the element.
% * ndf: Number of degrees of freedom (dof) on each face of the element.
% * xdof: Coordinates of the dof on the reference element. First row
% coincides with coordinate x_1 and second row with coordinate x_2.
% * xdof1d: Coordinates of the dof on the reference face (1d: edge).
% * Mref, Dxref, Dyref: Mass, stiffness matrices on the reference element.
% * invMref: Inverse of the mass matrix on the reference element.
% * M1d: Mass matrix on a reference segment [0, 1].
% * dof_face: Table of Nb of faces rows and ndf columns. For each face,
% it gives the indices of the dof on this face in the default orientation
% order.

switch parameter.pOrder
    
    case 0
        % P_0 case. p(x) = c
        % o
        % |`
        % | 1`
        % o----o
        nde = 1; ndf = 1;
        Pk_ref = struct('nde', nde, ...
            'ndf', ndf, ...
            'xdof', [1/3, 1/3]', ...
            'xdof1d', 1/3, ...
            'Mref', 1/2, ...
            'invMref', 2, ...
            'Dxref', 0, ...
            'Dyref', 0, ...
            'phie', @(xr, yr)1, ... % is not used in this HDG
            'M1d', 1, ...
            'invM1d', 1,...
            'dof_face', [1; 1; 1], ...
            'phi_face', @(s) 1 ...
            );
    case 1
        % P_1 case. p(x) = ax + by + c
        % 3
        % |`
        % |  `
        % 1----2
        nde = 3; ndf = 2;
        Pk_ref = struct('nde', nde, ...
            'ndf', ndf, ...
            'xdof', [0, 0; 1, 0; 0, 1]', ...
            'xdof1d', [0, 1], ...
            'Mref', 1./24.*...
            [2, 1, 1; 1, 2, 1; 1, 1, 2], ...
            'invMref', 2.*...
            [9, -3, -3; -3, 9, -3; -3, -3, 9], ...
            'Dxref', 1./6.*...
            [-1, -1, -1; 1, 1, 1; 0, 0, 0], ...
            'Dyref', 1./6.*...
            [-1, -1, -1; 0, 0, 0; 1, 1, 1], ...
            'phie', @(xr, yr)[1-xr-yr; xr; yr], ...
            'M1d', 1./6.*...
            [2, 1; 1, 2], ...
            'invM1d',...
            [4 -2; -2, 4],...
            'dof_face', [1, 2; 2, 3; 3, 1], ...
            'phi_face', @(s)[1-s; s] ...
            );
    case 2
        % P_2 case. p(x,y) = ax^2 + by^2 + cxy + dx + ey + f
        % order = 2, what we use. fig 3-4 in the PhD thesis.
        % 3
        % |`
        % 6  5
        % |    `
        % 1---4--2
        nde = 6; ndf = 3;
        Pk_ref = struct('nde', nde, ...
            'ndf', ndf, ...
            'xdof', 1./2.*...
            [0, 0; 2, 0; 0, 2; 1, 0; 1, 1; 0, 1]', ...
            'xdof1d', 1./2.*[0, 1, 2], ...
            'Mref', 1./360.*...
            [6, -1, -1,  0, -4,  0;
            -1,  6, -1,  0,  0, -4;
            -1, -1,  6, -4,  0,  0;
             0,  0, -4, 32, 16, 16;
            -4,  0,  0, 16, 32, 16;
             0, -4,  0, 16, 16, 32], ...
            'invMref', 1./4*...
            [288,  48,  48, -12,  48, -12;
              48, 288,  48, -12, -12,  48;
              48,  48, 288,  48, -12, -12;
             -12, -12,  48,  78, -27, -27;
              48, -12, -12, -27,  78, -27;
             -12, 48, -12,  -27, -27,  78], ...
            'Dxref', 1./30.*...
            [-2, -1, 0,  3, -1,  1;
              1,  2, 0, -3, -1,  1;
              1, -1, 0,  0,  2, -2;
             -3,  3, 0,  0,  4, -4;
              1,  3, 0, -4,  8, -8;
             -3, -1, 0,  4,  8, -8]', ...
            'Dyref', 1./30.*...
            [-2, 0, -1,  1, -1,  3;
              1, 0, -1, -2,  2,  0;
              1, 0,  2,  1, -1, -3;
             -3, 0, -1, -8,  8,  4;
              1, 0,  3, -8,  8, -4;
             -3, 0, 3,  -4,  4, 0]', ...
            'phie', @(xr, yr)[1-3*xr-3*yr+2*xr.^2+4*xr.*yr+2*yr.^2;
                              -1*xr+2*xr.^2;
                              -1*yr+2*yr.^2;
                              4*xr-4*xr.^2-4*xr.*yr;
                              4*xr.*yr;
                              4*yr-4*xr.*yr-4*yr.^2], ...
            'M1d', 1./30.*...
            [4, 2, -1; 2, 16, 2; -1, 2, 4], ...
            'invM1d',...
            [9, -3./2., 3; -3./2., 9./4., -3./2.; 3, -3./2., 9],...
            'dof_face', [1, 4, 2; 2, 5, 3; 3, 6, 1], ...
            'phi_face', @(s)[2*(1-s).*(.5-s); 4*s.*(1-s); 2*s.*(s-.5)] ...
            );
    case 3
        % P_3 case.
        % 3
        % | `
        % 8  7
        % |   `
        % 9 10 6
        % |     `
        % 1-4--5-2
        nde = 10; ndf = 4;
        Pk_ref = struct('nde', nde, ...
            'ndf', ndf, ...
            'xdof', 1./3.*...
            [0, 0; 3, 0; 0, 3; 1, 0; 2, 0; 2, 1; 1, 2; 0, 2; 0, 1; 1, 1]', ...
            'xdof1d', 1./3.*[0, 1, 2, 3], ...
            'Mref', .5/6720.*...
            [76.,   11.,   11.,   18.,    0.,   27.,   27.,    0.,   18.,   36.;
             11.,   76.,   11.,    0.,   18.,   18.,    0.,   27.,   27.,   36.;
             11.,   11.,   76.,   27.,   27.,    0.,   18.,   18.,    0.,   36.;
             18.,    0.,   27.,  540., -189., -135.,  -54., -135.,  270.,  162.;
              0.,   18.,   27., -189.,  540.,  270., -135.,  -54., -135.,  162.;
             27.,   18.,    0., -135.,  270.,  540., -189., -135.,  -54.,  162.;
             27.,    0.,   18.,  -54., -135., -189.,  540.,  270., -135.,  162.;
              0.,   27.,   18., -135.,  -54., -135.,  270.,  540., -189.,  162.;
             18.,   27.,    0.,  270., -135.,  -54., -135., -189.,  540.,  162.;
             36.,   36.,   36.,  162.,  162.,  162.,  162.,  162.,  162., 1944.], ...
            'invMref', 2./729.*...
            [7290., -729., -729., -297.,   54., -729., -729.,   54., -297.,   54.;
             -729., 7290., -729.,   54., -297., -297.,   54., -729., -729.,   54.;
             -729., -729., 7290., -729., -729.,   54., -297., -297.,   54.,   54.;
             -297.,   54., -729., 1570.,  455.,  367.,  218.,  367., -385., -198.;
               54., -297., -729.,  455., 1570., -385.,  367.,  218.,  367., -198.;
             -729., -297.,   54.,  367., -385., 1570.,  455.,  367.,  218., -198.;
             -729.,   54., -297.,  218.,  367.,  455., 1570., -385.,  367., -198.;
               54., -729., -297.,  367.,  218.,  367., -385., 1570.,  455., -198.;
             -297., -729.,   54., -385.,  367.,  218.,  367.,  455., 1570., -198.;
               54.,   54.,   54., -198., -198., -198., -198., -198., -198.,  348.], ...
            'Dxref', 1./3360.*...
            [-128.,   38.,    0.,  207., -117.,  45.,   45.,  -45.,    9.,  -54.;
              -38.,  128.,    0.,  117., -207.,  -9.,   45.,  -45.,  -45.,   54.;
              -38.,   38.,    0.,    0.,    0., -72.,  198., -198.,   72.,    0.;
             -207., -117.,    0.,    0.,  324.,-162.,  -81.,   81., -324.,  486.;
              117.,  207.,    0., -324.,    0., 324.,  -81.,   81.,  162., -486.;
              -45.,  207.,    0.,  162., -324., 648., -243.,  243.,  162., -810.;
              -45., -117.,    0.,   81.,   81.,  81.,  648., -648.,  243., -324.;
              117.,   45.,    0.,  -81.,  -81.,-243.,  648., -648.,  -81.,  324.;
             -207.,   45.,    0.,  324., -162.,-162., -243.,  243., -648.,  810.;
               54.,  -54.,    0., -486.,  486., 810.,  324., -324., -810.,  0.]', ...
            'Dyref', 1./3360.*...
            [-128.,    0.,   38.,    9.,  -45.,   45.,   45., -117.,  207.,  -54.;
              -38.,    0.,   38.,   72., -198.,  198.,  -72.,    0.,    0.,    0.;
              -38.,    0.,  128.,  -45.,  -45.,   45.,   -9., -207.,  117.,   54.;
             -207.,    0.,   45., -648.,  243., -243., -162., -162.,  324.,  810.;
              117.,    0.,   45.,  -81., -648.,  648., -243.,  -81.,  -81.,  324.;
              -45.,    0., -117.,  243., -648.,  648.,   81.,   81.,   81., -324.;
              -45.,    0.,  207.,  162.,  243., -243.,  648., -324.,  162., -810.;
              117.,    0.,  207.,  162.,   81.,  -81.,  324.,    0., -324., -486.;
             -207.,    0., -117., -324.,   81.,  -81., -162.,  324.,    0.,  486.;
              54.,    0.,  -54., -810., -324.,  324.,  810.,  486., -486.,    0.]', ...
            'M1d', 1./1680.*...
            [128., 99., -36., 19.;
             99., 648., -81., -36;
            -36., -81., 648., 99.;
             19., -36., 99., 128.], ...
            'invM1d',...
            [    16,   -68/27,    32/27,     -4;
             -68/27, 2224/729,   44/729,  32/27;
              32/27,   44/729, 2224/729, -68/27;
                 -4,    32/27,   -68/27,     16],...
            'dof_face', [1, 4, 5, 2; 2, 6, 7, 3; 3, 8, 9, 1], ...
            'phi_face', @(s)[-9./2.*(s-1./3.).*(s-2./3.).*(s-1);
                             27./2.*s.*(s-2./3.).*(s-1);
                             -27./2.*s.*(s-1./3.).*(s-1);
                             9./2.*s.*(s-1./3.).*(s-2./3.);
                            ] ...
            );       
    otherwise
        error('This polynomial order is not implemented.')
        
end

return;
function totMsh = gTriMesh(parameter)
% Load mesh

if parameter.prob == 1
    load disk_06_initial
elseif parameter.prob == 2
    load disk_06_ref1
elseif parameter.prob == 3
    load disk_06_ref2
elseif parameter.prob == 4
    load disk_06_ref3
end
totMsh.nCoord = p';%test coord
totMsh.eNodes = t';
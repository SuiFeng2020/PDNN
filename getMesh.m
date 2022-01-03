function totMsh = getMesh(parameter)
% read mesh file from gmsh and convert it into new style
if parameter.prob == 1 %
    filename = 'cylinder_16.msh';
elseif parameter.prob == 2
    filename = 'cylinder_26_1.msh';
elseif parameter.prob == 3 
    filename = 'cylinder_26_2.msh'; 
elseif parameter.prob == 4 
    filename = 'cylinder_26_3.msh'; 
else
    error('error in getMesh. The problem is not defined!')
end
%%
fid = fopen(filename, 'r');
%
buf = fscanf(fid, '%s', 1); % not interested
buf = fscanf(fid, '%s %i %i', 3); %
buf = fscanf(fid, '%s', 1); %
buf = fscanf(fid, '%s', 1); %
%
nNodes = fscanf(fid, '%i', 1); % number of nodes
nodes = fscanf(fid, '%i %f %f %f', [4, nNodes]); % All the nodes
nodes = nodes(2:3,:)';
%
buf = fscanf(fid, '%s', 1); % not interested
buf = fscanf(fid, '%s', 1); %
%
nElem = fscanf(fid, '%i', 1); % number of elements including boundary edges
% elements, first three rows are nodes, 4th row -- phyical region
% 5th row -- belong to which subdomain, 6th row -- possible neighboring
% subdomain if nonzero, if there are several?
elem = zeros(6,nElem);
elemTag = zeros(100,1);
bEdge = [];
n = 0; % number of elements
%
for i = 1 : nElem
    elemInfo = fscanf(fid, '%i %i %i', [3, 1]);
    for ii = 1 : elemInfo(3)
        elemTag(ii) = fscanf(fid, '%i', 1);
    end
    % boundary edges, maybe not useful
    if elemInfo(2) == 1
        edgeNodes = fscanf(fid,'%i %i',[2,1]);
        bEdge = [bEdge edgeNodes];
    end
    % triangle elements
    if elemInfo(2) == 2
        n = n + 1;
        elemNodes = fscanf(fid,'%i %i %i', [3, 1]);
        elem(1:3,n) = elemNodes;
        elem(4,n) = elemTag(1);
        elem(5,n) = elemTag(4);
        if elemInfo(3) > 4
            elem(6,n) = elemTag(5);
        end
    end
end
%
elems = elem(:,1:n)';
%
totMsh.nCoord = nodes;  % coordinates of nodes
totMsh.eNodes = elems;  % total elements

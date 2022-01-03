function dofmat = getdof
%
global totMsh totEdg refTri;
% Number of elements.
nElem = size(totMsh.eNodes,1);
% Computation  of the length of each face (edge here).
lenFaces = sqrt(...
    (totMsh.nCoord(totEdg.fNodes(:,2),1)-totMsh.nCoord(totEdg.fNodes(:,1),1)).^2 +...
    (totMsh.nCoord(totEdg.fNodes(:,2),2)-totMsh.nCoord(totEdg.fNodes(:,1),2)).^2);
dofmat.hmin = min(lenFaces);
dofmat.hmax = max(lenFaces);
% Record the DOFs
dofmat.DOF = zeros(nElem, refTri.nde,2);
Xdof = zeros(nElem, refTri.nde);
Ydof = zeros(nElem, refTri.nde);
dofmat.absdetje = zeros(nElem, 1); 
dofmat.dtscale = dofmat.absdetje;
for k = 1:nElem 
    jac = [totMsh.nCoord(totMsh.eNodes(k,2),1)-totMsh.nCoord(totMsh.eNodes(k,1),1),...
           totMsh.nCoord(totMsh.eNodes(k,3),1)-totMsh.nCoord(totMsh.eNodes(k,1),1);...
           totMsh.nCoord(totMsh.eNodes(k,2),2)-totMsh.nCoord(totMsh.eNodes(k,1),2),...
           totMsh.nCoord(totMsh.eNodes(k,3),2)-totMsh.nCoord(totMsh.eNodes(k,1),2)];
    len1 = sqrt((totMsh.nCoord(totMsh.eNodes(k,2),1)-totMsh.nCoord(totMsh.eNodes(k,1),1))^2 + ...
           (totMsh.nCoord(totMsh.eNodes(k,2),2)-totMsh.nCoord(totMsh.eNodes(k,1),2))^2);
    len2 = sqrt((totMsh.nCoord(totMsh.eNodes(k,3),1)-totMsh.nCoord(totMsh.eNodes(k,2),1))^2 + ...
           (totMsh.nCoord(totMsh.eNodes(k,3),2)-totMsh.nCoord(totMsh.eNodes(k,2),2))^2);
    len3 = sqrt((totMsh.nCoord(totMsh.eNodes(k,1),1)-totMsh.nCoord(totMsh.eNodes(k,3),1))^2 + ...
           (totMsh.nCoord(totMsh.eNodes(k,1),2)-totMsh.nCoord(totMsh.eNodes(k,3),2))^2);
    %p=(a+b+c)/2 S=sqrt[p(p-a)(p-b)(p-c)] 
    sper = (len1 + len2 + len3)/2.0; 
    Area = sqrt(sper.*(sper-len1).*(sper-len2).*(sper-len3));
%     if totMsh.eNodes(k,4) == Metal
%         localc = 1/sqrt(1*eps);
%     elseif totMsh.eNodes(k,4) == FreeSpac
%         localc = 1/sqrt(1*1);
%     end
%     % Compute scale using radius of inscribed circle
%     dofmat.dtscale(k) = Area./(localc*sper*2);
    detj = det(jac);
    dofmat.absdetje(k) = abs(detj); % twice of the area of the triangle
    Tmp = jac*refTri.xdof+...
        kron(ones(1, refTri.nde), totMsh.nCoord(totMsh.eNodes(k, 1),:)');
    Xdof(k, :) = Tmp(1, :);%dof 1 2 3 4 5 6 but the edge 1 4 2;2 5 3;3 6 1;
    Ydof(k, :) = Tmp(2, :);
end      
dofmat.DOF(:,:,1) = Xdof;
dofmat.DOF(:,:,2) = Ydof;
[mapM, mapP, dofmat.vmapM, dofmat.vmapP, dofmat.vmapB, dofmat.mapB,dofmat.EToE] = ...
                    buildMaps(totMsh.eNodes,totMsh.nCoord,refTri,dofmat.DOF);
end  
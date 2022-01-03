function ADGTDinvMeeps = getinvMe(eps)
%
global totMsh refTri;
% Number of elements.
% test for the computational domain property
ind = find(totMsh.eNodes(:,4)==1); n1 = length(ind);
ind = find(totMsh.eNodes(:,4)==2); n2 = length(ind);
if n1 < n2
    Metal = 1; FreeSpac = 2;
else
    Metal = 2; FreeSpac = 1;
end
nElem = size(totMsh.eNodes,1);
Ndof = refTri.nde*nElem;
nodeids = reshape(1:nElem*refTri.nde,refTri.nde,nElem);
absdetje = zeros(nElem, 1);
%the mass maxtrix
ADGTDinvMeeps = sparse(nElem*refTri.nde,nElem*refTri.nde);
A.invMeeps = zeros(nElem,refTri.nde,refTri.nde);
indi = zeros(nElem,refTri.nde,refTri.nde);
indj = zeros(nElem,refTri.nde,refTri.nde);
for k = 1:nElem  
    jac = [totMsh.nCoord(totMsh.eNodes(k,2),1)-totMsh.nCoord(totMsh.eNodes(k,1),1), ...
           totMsh.nCoord(totMsh.eNodes(k,3),1)-totMsh.nCoord(totMsh.eNodes(k,1),1);...
           totMsh.nCoord(totMsh.eNodes(k,2),2)-totMsh.nCoord(totMsh.eNodes(k,1),2), ...
           totMsh.nCoord(totMsh.eNodes(k,3),2)-totMsh.nCoord(totMsh.eNodes(k,1),2)];
    detj = det(jac);
    absdetje(k) = abs(detj); % twice of the area of the triangle
    inde = nodeids(:,k);
    % identify the location of the triangle
    if totMsh.eNodes(k,4) == FreeSpac
        ieps = 1;
    elseif totMsh.eNodes(k,4) == Metal
        ieps = eps;
    end
    invMeeeps    = (1/(ieps*absdetje(k)))*refTri.invMref;
    A.invMeeps(k,:,:) = invMeeeps;
    indi(k,:,:) = kron(ones(1,refTri.nde),inde);
    indj(k,:,:) = squeeze(indi(k,:,:))';
end   
%
% assembing the golbal matrix
ADGTDinvMeeps = sparse(reshape(indi, numel(indi),1),reshape(indj, numel(indj),1),...
    reshape(A.invMeeps, numel(A.invMeeps),1), Ndof, Ndof);
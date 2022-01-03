function totEdg = getFace
% [fNodes, fElem, eFaces] = getFace(eNodes)
%
% Find out all the faces and the associated elements
% Inputs:    eNodes -- matrix(nT,3) containing nodes of the tetras
%
% Outputs:   fNodes -- matrix(nF,2) containing nodes of the faces
%            fElem  -- matrix(nF,2) elements associated with each face
%                      negativeness means outward normal of the element
%                      has opposite sign with the face
%            eFaces -- matrix(nElem,3) elements represented by faces

% Author: PLUM, STEPHANE
% Date  : 2014-07-30

% Note that we have some redutant information, but it's not a problem
global totMsh;
%% memory allocation and initialization
nbElem = size(totMsh.eNodes,1);
totEdg.fNodes = zeros(3*nbElem,2);
totEdg.fElem = zeros(4*nbElem,2);
totEdg.eFaces = zeros(nbElem,3);%
nF = 0; 
%% Loop over all the triangles
for k = 1 : nbElem
    % Nodes and faces are counter-clockwise orienated
    % find faces from faceK(3,2)
    faceK = [totMsh.eNodes(k,[1,2]);totMsh.eNodes(k,[2,3]);totMsh.eNodes(k,[3,1])];
    for j = 1 : 3 % each triangle has 3 faces
        % the use of function 'find' makes it slow
        indF = find(totEdg.fNodes(1:nF,:) == faceK(j,1));
        if (isempty(indF)) % find a new face, because this node has not appeared
            nF = nF + 1;
            totEdg.fNodes(nF,:) = faceK(j,:);
            totEdg.fElem(nF,1) = k;
            totEdg.eFaces(k,j) = nF;
%            eFacePerm(k,j) = 1;
        else
            indF = mod(indF-1,nF)+1;
            lenIndF = length(indF);
            flg = 0;
            for l = 1 : lenIndF
                if (isempty(setdiff(faceK(j,:), totEdg.fNodes(indF(l),:)))) 
                    % the face has appeared
                    % the second element associated with this face is assigned with a minus sign
                    totEdg.fElem(indF(l),2) = -k; 
                    tmpNodes = totEdg.fNodes(indF(l),:);
                    % minus sign states the orientation
                    totEdg.eFaces(k,j) = -indF(l);
                    flg = 1;
                    break;
                end
            end
            if (flg == 0) % find a new face
                nF = nF + 1;
                totEdg.fNodes(nF,:) = faceK(j,:);
                totEdg.fElem(nF,1) = k;
                totEdg.eFaces(k,j) = nF;
            end
        end
    end
end
% trim the results
totEdg.fNodes = totEdg.fNodes(1:nF,:);
totEdg.fElem = totEdg.fElem(1:nF, :);
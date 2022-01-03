function [errE, errH] = getErr(poddgsolution,dgsolution)
%
global ADGTD;
% calculation of the error
poddgH = poddgsolution(:,1);
poddgE = poddgsolution(:,2);

dgH = dgsolution(:,1);
dgE = dgsolution(:,2);
errH = zeros(1,1);
errE = zeros(1,1);
% for ii = 1:2
%    errH(ii,1) = (dgH(:,ii) - poddgH(:,ii))'*ADGTD.Me...
%                                            *(dgH(:,ii) - poddgH(:,ii));
% end
% 
% for ii = 1:1
% errE(ii,1) = (dgE(:,ii) - poddgE(:,ii))'*ADGTD.Me...
%                                           *(dgE(:,ii) - poddgE(:,ii));
% end
% euclidean error
for ii = 1:1
   errH(ii,1) = (dgH(:,ii) - poddgH(:,ii))'...
                                           *(dgH(:,ii) - poddgH(:,ii));
end

for ii = 1:1
errE(ii,1) = (dgE(:,ii) - poddgE(:,ii))'...
                                          *(dgE(:,ii) - poddgE(:,ii));
end
errE = sum(sqrt(real(errE)));
errH = sum(sqrt(real(errH)));
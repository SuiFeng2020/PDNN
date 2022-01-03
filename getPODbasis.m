function POD = getPODbasis(info,Snapshots)
% the calculation of the POD mehod
% Method 1 : SVD of A (snapshots matrix)
% Method 2 : eig of AA^T
% Method 3 : eig of A^TA
%% Calculating the eigs based on the method 3
[Vect.Hx, Val.Hx] = eig(Snapshots.Hxe'*Snapshots.Hxe);
[Vect.Hy, Val.Hy] = eig(Snapshots.Hye'*Snapshots.Hye);
[Vect.Ez, Val.Ez] = eig(Snapshots.Eze'*Snapshots.Eze);

Vect.Hx = Vect.Hx(:, end:-1:1);
Val.Hx = abs(Val.Hx(end:-1:1,end:-1:1));
Vect.Hy = Vect.Hy(:, end:-1:1);
Val.Hy = abs(Val.Hy(end:-1:1,end:-1:1));
Vect.Ez = Vect.Ez(:, end:-1:1);
Val.Ez = abs(Val.Ez(end:-1:1,end:-1:1));

% display the variation of the sigular values
% figure(1)
% semilogy(1:size(diag(Val.Hx)), diag(sqrt(Val.Hx)),'-ro',...
%          'LineWidth',1)
% hold on
% semilogy(1:size(diag(Val.Hy)), diag(sqrt(Val.Hy)),'-go',...
%          'LineWidth',1)
% 
% semilogy(1:size(diag(Val.Ez)), diag(sqrt(Val.Ez)),'k<-',...
%          'LineWidth',1)
% grid on
% legend('Hx','Hy','Ez')

% Computing the dimention of ROM based on a energy mehod
val.Hx = diag(Val.Hx);
val.Hy = diag(Val.Hy);
val.Ez = diag(Val.Ez);

% 
elli = 1;
while (sum(val.Hx(1:elli))/sum(val.Hx) < info)
    elli = elli + 1;
end
POD.Dimen.Hx = elli;

elli = 1;
while (sum(val.Hy(1:elli))/sum(val.Hy) < info)
    elli = elli + 1;
end
POD.Dimen.Hy = elli;

elli = 1;
while (sum(val.Ez(1:elli))/sum(val.Ez) < info)
    elli = elli + 1;
end
POD.Dimen.Ez = elli;

% for ii = 1:40
%     Cumulative_energy.Hx(ii) = sum(val.Hx(1:ii))/sum(val.Hx);
% end
% 
% for ii = 1:40
%     Cumulative_energy.Hy(ii) = sum(val.Hy(1:ii))/sum(val.Hy);
% end
% 
% for ii = 1:40
%     Cumulative_energy.Ez(ii) = sum(val.Ez(1:ii))/sum(val.Ez);
% end

% figure(2)
% h = gca;
% set(h,'FontSize',20);
% plot(1:40, Cumulative_energy.Hx,'-ro',...
%          'LineWidth',2)
% hold on
% plot(1:40, Cumulative_energy.Hy,'-k+',...
%          'LineWidth',2)
% 
% plot(1:40, Cumulative_energy.Ez,'-bx',...
%          'LineWidth',2)
% xlabel('Number of Modes ');
% ylabel('Cumulative Energy of the Eigenvalues');
% legend('Hx','Hy','Ez')
% grid on
% set(gcf,'PaperPositionmode','auto')

% POD.Dimen.Hx = sum(sum(val.Hx~=0));
% POD.Dimen.Hy = sum(sum(val.Hy~=0));
% POD.Dimen.Ez = sum(sum(val.Ez~=0));

% Getting the POD basis
% for ii = 1:POD.Dimen.Hx
%     POD.Basis.Hx(:,ii) = Snapshots.Hxe*Vect.Hx(:,ii)...
%                         /sqrt(Val.Hx(ii,ii));
% end
% for ii = 1:POD.Dimen.Hy
%     POD.Basis.Hy(:,ii) = Snapshots.Hye*Vect.Hy(:,ii)...
%                         /sqrt(Val.Hy(ii,ii));                  
% end
% 
% for ii = 1:POD.Dimen.Ez
%     POD.Basis.Ez(:,ii) = Snapshots.Eze*Vect.Ez(:,ii)...
%                         /sqrt(Val.Ez(ii,ii));
% end
Ndof = size(Snapshots.Hxe,1);
POD.Basis.Hx = (Snapshots.Hxe*Vect.Hx(:,1:POD.Dimen.Hx))...
                                 ./(ones(Ndof,1)*sqrt(val.Hx(1:POD.Dimen.Hx,1))');
POD.Basis.Hy = (Snapshots.Hye*Vect.Hy(:,1:POD.Dimen.Hy))...
                                ./(ones(Ndof,1)*sqrt(val.Hy(1:POD.Dimen.Hy,1))');                  
POD.Basis.Ez = (Snapshots.Eze*Vect.Ez(:,1:POD.Dimen.Ez))...
                               ./(ones(Ndof,1)*sqrt(val.Ez(1:POD.Dimen.Ez,1))');
%% Calculating the eigs based on the method 1
% comparison with method 3 
% [Vect.Hx, Val.Hx] = svd(Snapshots.Hxe,0);
% [Vect.Hy, Val.Hy] = svd(Snapshots.Hye,0);
% [Vect.Ez, Val.Ez] = svd(Snapshots.Eze,0);
% 
% val.Hx = diag(Val.Hx);
% val.Hy = diag(Val.Hy);
% val.Ez = diag(Val.Ez);
% 
% % 
% elli = 1;
% while (sum((val.Hx(1:elli)).^2)/sum((val.Hx).^2) < info)
%     elli = elli + 1;
% end
% POD.Dimen.Hx = elli;
% 
% elli = 1;
% while (sum((val.Hy(1:elli)).^2)/sum((val.Hy).^2) < info)
%     elli = elli + 1;
% end
% POD.Dimen.Hy = elli;
% 
% elli = 1;
% while (sum((val.Ez(1:elli)).^2)/sum((val.Ez).^2) < info)
%     elli = elli + 1;
% end
% POD.Dimen.Ez = elli;
% % Getting the POD basis
% POD.Basis.Hx = Vect.Hx(:,1:POD.Dimen.Hx);
% POD.Basis.Hy = Vect.Hy(:,1:POD.Dimen.Hy);                  
% POD.Basis.Ez = Vect.Ez(:,1:POD.Dimen.Ez);
end

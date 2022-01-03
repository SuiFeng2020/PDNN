function showresultsANN
% 1 5 10 ( 2 4 6 8)
clc;
global fphie fpsie fphihy fpsihy fphihx fpsihx;
global fphiephieltol fpsiephieltol fphihyphihyltol fpsihyphihyltol ...
       fphihxphihxltol fpsihxphihxltol;
global traindata;
c = [0 0 0;1 0 0;0.5 0.16 0.16;0 0 1;0 1 0;];
xtest = mapminmax(traindata.time);
ptest = mapminmax(traindata.parameter);
% xtest = traindata.timeNor;
% ptest = traindata.parameterNor;
l1 = 1; l2 = 5; l3 = 10;
kk = 2:2:10;
figure();
l = l1;
cii = 1;
hh = subplot(3,2,1);
ppp = get(hh,'pos');
ppp(3) = ppp(3) + 0.04;
ppp(4) = ppp(4) - 0.05;
set(hh,'pos',ppp)
for k = kk
    [pred,~,ci] = sim(fphie.(['l_' num2str(l) 'qle' num2str(k)]),xtest);
    h = gca;
    set(h,'FontSize',15);
    phieltollk = fphiephieltol.(['l_' num2str(l)  'qle' num2str(k)]);
    h1 = plot(traindata.time',phieltollk,'color',c(cii,:),'LineWidth',1.5);
    hold on
    h2 = plot(xtest(1:3:end),pred(1:3:end),'color',c(cii,:),'LineWidth',1.5);
    set(h1,'LineStyle','-')
    set(h2,'LineStyle','--','Marker','+','MarkerSize',5)
    hold on
    cii = cii + 1;
end
xlim([min(traindata.time),max(traindata.time)])
xlabel('Time (m)')
grid on
title('The 1th coefficent - time modes')
%
cii = 1;
hh = subplot(3,2,2);
ppp = get(hh,'pos');
ppp(3) = ppp(3) + 0.04;
ppp(4) = ppp(4) - 0.05;
set(hh,'pos',ppp)
for k = kk
    [pred,~,ci] = sim(fpsie.(['l_' num2str(l) 'qle' num2str(k)]),ptest);
    h = gca;
    set(h,'FontSize',15);
    psieltollk = fpsiephieltol.(['l_' num2str(l)  'qle' num2str(k)]);
    h1 = plot(traindata.parameter(1:1:end)',psieltollk(1:1:end),'color',c(cii,:),'LineWidth',1.5);
    hold on
    h2 = plot(ptest(1:1:end),pred(1:1:end),'color',c(cii,:),'LineWidth',1.5);
    set(h1,'LineStyle','-')
    set(h2,'LineStyle','--','Marker','+','MarkerSize',5)
    hold on
    cii = cii + 1;
end
xlim([min(traindata.parameter),max(traindata.parameter)])
xlabel('$$\varepsilon$$','Interpreter','Latex')
grid on
title('The 1th coefficent - parameter modes')
%%
l = l2;
cii = 1;
hh = subplot(3,2,3);
ppp = get(hh,'pos');
ppp(3) = ppp(3) + 0.04;
ppp(4) = ppp(4) - 0.05;
set(hh,'pos',ppp)
for k = kk
    [pred,~,ci] = sim(fphie.(['l_' num2str(l) 'qle' num2str(k)]),xtest);
    h = gca;
    set(h,'FontSize',15);
    phieltollk = fphiephieltol.(['l_' num2str(l)  'qle' num2str(k)]);
    h1 = plot(traindata.time',phieltollk,'color',c(cii,:),'LineWidth',1.5);
    hold on
    h2 = plot(xtest(1:3:end),pred(1:3:end),'color',c(cii,:),'LineWidth',1.5);
    set(h1,'LineStyle','-')
    set(h2,'LineStyle','--','Marker','+','MarkerSize',5)
    hold on
    cii = cii + 1;
end
xlim([min(traindata.time),max(traindata.time)])
xlabel('Time (m)')
grid on
title('The 5th coefficent - time modes')
%
cii = 1;
hh = subplot(3,2,4);
ppp = get(hh,'pos');
ppp(3) = ppp(3) + 0.04;
ppp(4) = ppp(4) - 0.05;
set(hh,'pos',ppp)
for k = 2:2:10
    [pred,~,ci] = sim(fpsie.(['l_' num2str(l) 'qle' num2str(k)]),ptest);
    h = gca;
    set(h,'FontSize',15);
    psieltollk = fpsiephieltol.(['l_' num2str(l)  'qle' num2str(k)]);
    h1 = plot(traindata.parameter',psieltollk,'color',c(cii,:),'LineWidth',1.5);
    hold on
    h2 = plot(ptest(1:1:end),pred(1:1:end),'color',c(cii,:),'LineWidth',1.5);
    set(h1,'LineStyle','-')
    set(h2,'LineStyle','--','Marker','+','MarkerSize',5)
    hold on
    cii = cii + 1;
end
xlim([min(traindata.parameter),max(traindata.parameter)])
xlabel('$$\varepsilon$$','Interpreter','Latex')
grid on
title('The 5th coefficent - parameter modes')
%%
%
l = l3;
cii = 1;
hh = subplot(3,2,5);
ppp = get(hh,'pos');
ppp(3) = ppp(3) + 0.04;
ppp(4) = ppp(4) - 0.05;
set(hh,'pos',ppp)
for k = kk
    [pred,~,ci] = sim(fphie.(['l_' num2str(l) 'qle' num2str(k)]),xtest);
    h = gca;
    set(h,'FontSize',15);
    phieltollk = fphiephieltol.(['l_' num2str(l)  'qle' num2str(k)]); 
    h1 = plot(traindata.time',phieltollk,'color',c(cii,:),'LineWidth',1.5);
    hold on
    h2 = plot(xtest(1:3:end),pred(1:3:end),'color',c(cii,:),'LineWidth',1.5);
    set(h1,'LineStyle','-')
    set(h2,'LineStyle','--','Marker','+','MarkerSize',5)
    hold on
    cii = cii + 1;
end
xlim([min(traindata.time),max(traindata.time)])
xlabel('Time (m)')
grid on
title('The 10th coefficent - time modes')
%
cii = 1;
hh = subplot(3,2,6);
ppp = get(hh,'pos');
ppp(3) = ppp(3) + 0.04;
ppp(4) = ppp(4) - 0.05;
set(hh,'pos',ppp)
for k =kk
    [pred,~,ci] = sim(fpsie.(['l_' num2str(l) 'qle' num2str(k)]),ptest);
    h = gca;
    set(h,'FontSize',15);
    psieltollk = fpsiephieltol.(['l_' num2str(l)  'qle' num2str(k)]);
    h1 = plot(traindata.parameter',psieltollk,'color',c(cii,:),'LineWidth',1.5);
    hold on
    h2 = plot(ptest(1:1:end),pred(1:1:end),'color',c(cii,:),'LineWidth',1.5);
    set(h1,'LineStyle','-')
    set(h2,'LineStyle','--','Marker','+','MarkerSize',5)
    hold on
    cii = cii + 1;
end
xlim([min(traindata.parameter),max(traindata.parameter)])
xlabel('$$\varepsilon$$','Interpreter','Latex')
grid on
title('The 10th coefficent - parameter modes')
%% Hy
figure();
l = l1;
cii = 1;
hh = subplot(3,2,1);
ppp = get(hh,'pos');
ppp(3) = ppp(3) + 0.04;
ppp(4) = ppp(4) - 0.05;
set(hh,'pos',ppp)
for k = kk
    [pred,~,ci] = sim(fphihy.(['l_' num2str(l) 'qlhy' num2str(k)]),xtest);
    h = gca;
    set(h,'FontSize',15);
    phihyltollk = fphihyphihyltol.(['l_' num2str(l)  'qlhy' num2str(k)]);
    h1 = plot(traindata.time',phihyltollk,'color',c(cii,:),'LineWidth',1.5);
    hold on
    h2 = plot(xtest(1:3:end),pred(1:3:end),'color',c(cii,:),'LineWidth',1.5);
    set(h1,'LineStyle','-')
    set(h2,'LineStyle','--','Marker','+','MarkerSize',5)
    hold on
    cii = cii + 1;
end
xlim([min(traindata.time),max(traindata.time)])
xlabel('Time (m)')
grid on
title('The 1th coefficent - time modes')
%
cii = 1;
hh = subplot(3,2,2);
ppp = get(hh,'pos');
ppp(3) = ppp(3) + 0.04;
ppp(4) = ppp(4) - 0.05;
set(hh,'pos',ppp)
for k = kk
    [pred,~,ci] = sim(fpsihy.(['l_' num2str(l) 'qlhy' num2str(k)]),ptest);
    h = gca;
    set(h,'FontSize',15);
    psihyltollk = fpsihyphihyltol.(['l_' num2str(l)  'qlhy' num2str(k)]);
    h1 = plot(traindata.parameter(1:1:end)',psihyltollk(1:1:end),'color',c(cii,:),'LineWidth',1.5);
    hold on
    h2 = plot(ptest(1:1:end),pred(1:1:end),'color',c(cii,:),'LineWidth',1.5);
    set(h1,'LineStyle','-')
    set(h2,'LineStyle','--','Marker','+','MarkerSize',5)
    hold on
    cii = cii + 1;
end
xlim([min(traindata.parameter),max(traindata.parameter)])
xlabel('$$\varepsilon$$','Interpreter','Latex')
grid on
title('The 1th coefficent - parameter modes')
%%
l = l2;
cii = 1;
hh = subplot(3,2,3);
ppp = get(hh,'pos');
ppp(3) = ppp(3) + 0.04;
ppp(4) = ppp(4) - 0.05;
set(hh,'pos',ppp)
for k = kk
    [pred,~,ci] = sim(fphihy.(['l_' num2str(l) 'qlhy' num2str(k)]),xtest);
    h = gca;
    set(h,'FontSize',15);
    phihyltollk = fphihyphihyltol.(['l_' num2str(l)  'qlhy' num2str(k)]);
    h1 = plot(traindata.time',phihyltollk,'color',c(cii,:),'LineWidth',1.5);
    hold on
    h2 = plot(xtest(1:3:end),pred(1:3:end),'color',c(cii,:),'LineWidth',1.5);
    set(h1,'LineStyle','-')
    set(h2,'LineStyle','--','Marker','+','MarkerSize',5)
    hold on
    cii = cii + 1;
end
xlim([min(traindata.time),max(traindata.time)])
xlabel('Time (m)')
grid on
title('The 5th coefficent - time modes')
%
cii = 1;
hh = subplot(3,2,4);
ppp = get(hh,'pos');
ppp(3) = ppp(3) + 0.04;
ppp(4) = ppp(4) - 0.05;
set(hh,'pos',ppp)
for k = 2:2:10
    [pred,~,ci] = sim(fpsihy.(['l_' num2str(l) 'qlhy' num2str(k)]),ptest);
    h = gca;
    set(h,'FontSize',15);
    psihyltollk = fpsihyphihyltol.(['l_' num2str(l)  'qlhy' num2str(k)]);
    h1 = plot(traindata.parameter',psihyltollk,'color',c(cii,:),'LineWidth',1.5);
    hold on
    h2 = plot(ptest(1:1:end),pred(1:1:end),'color',c(cii,:),'LineWidth',1.5);
    set(h1,'LineStyle','-')
    set(h2,'LineStyle','--','Marker','+','MarkerSize',5)
    hold on
    cii = cii + 1;
end
xlim([min(traindata.parameter),max(traindata.parameter)])
xlabel('$$\varepsilon$$','Interpreter','Latex')
grid on
title('The 5th coefficent - parameter modes')
%%
%
l = l3;
cii = 1;
hh = subplot(3,2,5);
ppp = get(hh,'pos');
ppp(3) = ppp(3) + 0.04;
ppp(4) = ppp(4) - 0.05;
set(hh,'pos',ppp)
for k = kk
    [pred,~,ci] = sim(fphihy.(['l_' num2str(l) 'qlhy' num2str(k)]),xtest);
    h = gca;
    set(h,'FontSize',15);
    phihyltollk = fphihyphihyltol.(['l_' num2str(l)  'qlhy' num2str(k)]); 
    h1 = plot(traindata.time',phihyltollk,'color',c(cii,:),'LineWidth',1.5);
    hold on
    h2 = plot(xtest(1:3:end),pred(1:3:end),'color',c(cii,:),'LineWidth',1.5);
    set(h1,'LineStyle','-')
    set(h2,'LineStyle','--','Marker','+','MarkerSize',5)
    hold on
    cii = cii + 1;
end
xlim([min(traindata.time),max(traindata.time)])
xlabel('Time (m)')
grid on
title('The 10th coefficent - time modes')
%
cii = 1;
hh = subplot(3,2,6);
ppp = get(hh,'pos');
ppp(3) = ppp(3) + 0.04;
ppp(4) = ppp(4) - 0.05;
set(hh,'pos',ppp)
for k =kk
    [pred,~,ci] = sim(fpsihy.(['l_' num2str(l) 'qlhy' num2str(k)]),ptest);
    h = gca;
    set(h,'FontSize',15);
    psihyltollk = fpsihyphihyltol.(['l_' num2str(l)  'qlhy' num2str(k)]);
    h1 = plot(traindata.parameter',psihyltollk,'color',c(cii,:),'LineWidth',1.5);
    hold on
    h2 = plot(ptest(1:1:end),pred(1:1:end),'color',c(cii,:),'LineWidth',1.5);
    set(h1,'LineStyle','-')
    set(h2,'LineStyle','--','Marker','+','MarkerSize',5)
    hold on
    cii = cii + 1;
end
xlim([min(traindata.parameter),max(traindata.parameter)])
xlabel('$$\varepsilon$$','Interpreter','Latex')
grid on
title('The 10th coefficent - parameter modes')
% print('-depsc','cylindergpr1510thtimemodes')
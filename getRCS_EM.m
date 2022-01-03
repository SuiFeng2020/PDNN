function RCS = getRCS_EM(dgsolution,podsolution,DOF,r_sca,r_imag,parameter)
% calculating RCS based on EM obtained the POD-DGTD and DGTD methods
% getRCS_EM(out.dgtdsolution,out.dgtdsolution,out.dof,0.6,1.4,parameter)
%%
% r_imag: the radius of the imaginary closed surface
% r_obser: the radius of the observation locations
st_sca = 0.005*pi;
theta_sca = 0:st_sca:2*pi-st_sca;
x_sca = r_sca*cos(theta_sca);
y_sca = r_sca*sin(theta_sca);
%
st_imag = 0.005*pi;
theta_imag = 0:st_imag:2*pi-st_imag;
x_imag = r_imag*cos(theta_imag);
y_imag = r_imag*sin(theta_imag);
nx_imag = cos(theta_imag);
ny_imag = sin(theta_imag);
len = 2*r_imag*sin(st_imag/2);
%
st_obser = 0.005*pi;
theta_obser = -pi:st_obser:pi;
rx_obser = cos(theta_obser);
ry_obser = sin(theta_obser);
%
omega = 2*pi*parameter.freq/parameter.c0; % normalied
k = sqrt(1*1)*omega/(parameter.c); % the normalized speed of c = 1;
%% Incident fields E^inc and H^inc at imaginary surface
tt = 0:parameter.dt:parameter.tmax*parameter.c0;
% Ezinc = @(x,y,t)(cos(omega*t - k*x));
% Hyinc = @(x,y,t)(-cos(omega*t - k*x));
Hyinc = @(x,y,t)(-exp((1i)*omega*t - (1i)*k*x));
Ezinc = @(x,y,t)(exp((1i)*omega*t - (1i)*k*x));
Ezinc_timesca = zeros(length(x_sca),length(tt));
Ezinc_imagtime = zeros(length(x_imag),length(tt));
Hyinc_imagtime = zeros(length(x_imag),length(tt));
for i = 1:length(tt)
    Ezinc_timesca(:,i) = Ezinc(x_sca,y_sca,tt(i));
    Ezinc_imagtime(:,i) = Ezinc(x_imag,y_imag,tt(i));
    Hyinc_imagtime(:,i) = Hyinc(x_imag,y_imag,tt(i));
end
T = parameter.dt;
Fs = 1/T;
L = length(tt);
% NFFT = 2^nextpow2(L);
Fx = 0:Fs/(L-1):Fs;
Fn = parameter.freq/parameter.c0;
fftf = fft(Ezinc_timesca');
Ezinc_scafreq = interp1(Fx',fftf,Fn)/(L);
absEzincsca = sqrt(sum(Ezinc_scafreq.*conj(Ezinc_scafreq)));
fftf = fft(Ezinc_imagtime');
Ezinc_imagfreq = interp1(Fx',fftf,Fn)/(L);
fftf = fft(Hyinc_imagtime');
Hyinc_imagfreq = interp1(Fx',fftf,Fn)/(L);
%% Scatter fields E^sca and H^sca at imaginary surface
xDod = DOF(:,:,1)';
yDod = DOF(:,:,2)';
[m,n] = size(xDod);nNod  = m*n;
xreshape = reshape(xDod,nNod,1);
yreshape = reshape(yDod,nNod,1);
%
Hxe.dg = dgsolution(:,1);
Hxe.pod = podsolution(:,1);
Hye.dg = dgsolution(:,2);
Hye.pod = podsolution(:,2);
Eze.dg = dgsolution(:,3);
Eze.pod = podsolution(:,3);
%
fHxdg = scatteredInterpolant(xreshape,yreshape,Hxe.dg);
fHxpod = scatteredInterpolant(xreshape,yreshape,Hxe.pod);
fHydg = scatteredInterpolant(xreshape,yreshape,Hye.dg);
fHypod = scatteredInterpolant(xreshape,yreshape,Hye.pod);
fEzdg = scatteredInterpolant(xreshape,yreshape,Eze.dg);
fEzpod = scatteredInterpolant(xreshape,yreshape,Eze.pod);
%
Hxsca_imagfreq.dg = fHxdg(x_imag,y_imag) - 0; 
Hxsca_imagfreq.pod = fHxpod(x_imag,y_imag) - 0;
Hysca_imagfreq.dg = fHydg(x_imag,y_imag) - Hyinc_imagfreq; 
Hysca_imagfreq.pod = fHypod(x_imag,y_imag) - Hyinc_imagfreq;
Ezsca_imagfreq.dg = fEzdg(x_imag,y_imag) - Ezinc_imagfreq; 
Ezsca_imagfreq.pod = fEzpod(x_imag,y_imag) - Ezinc_imagfreq;
%% Radar cross section (RCS)
RCS.dg = zeros(1,length(theta_obser));
RCS.pod = zeros(1,length(theta_obser));
for ii = 1:length(theta_obser) % for iith observation location
    for jj = 1:length(theta_imag) % summation
      sumiijj.dg(jj) = ...
          (parameter.c*1*(nx_imag(jj)*Hysca_imagfreq.dg(jj) ...
                     - ny_imag(jj)*Hxsca_imagfreq.dg(jj))...
        + (rx_obser(ii)*nx_imag(jj)*Ezsca_imagfreq.dg(jj)...
          + ry_obser(ii)*ny_imag(jj)*Ezsca_imagfreq.dg(jj)))...
        *exp(-(1i)*k*(rx_obser(ii)*nx_imag(jj) ...
                     + ry_obser(ii)*ny_imag(jj))*r_imag);
      sumiijj.pod(jj) = ...
          (parameter.c*1*(nx_imag(jj)*Hysca_imagfreq.pod(jj) ...
                     - ny_imag(jj)*Hxsca_imagfreq.pod(jj))...
        + (rx_obser(ii)*nx_imag(jj)*Ezsca_imagfreq.pod(jj)...
          + ry_obser(ii)*ny_imag(jj)*Ezsca_imagfreq.pod(jj)))...
        *exp(-(1i)*k*(rx_obser(ii)*nx_imag(jj) ...
                     + ry_obser(ii)*ny_imag(jj))*r_imag);            
    end
%     I.dg(ii) = sum(sumiijj.dg*len)/absEzincsca;
%     I.pod(ii) = sum(sumiijj.pod*len)/absEzincsca;
    I.dg(ii) = sum(sumiijj.dg*len);
    I.pod(ii) = sum(sumiijj.pod*len);
end
% RCS.dg = (k^2/(4*pi))*(abs(I.dg).^2);
% RCS.pod = (k^2/(4*pi))*(abs(I.pod).^2);
RCS.dg = (k/(4))*(abs(I.dg).^2);
RCS.pod = (k/(4))*(abs(I.pod).^2);
RCS.theta = theta_obser + pi;
% plot((180/pi)*(theta_obser+pi),10*log10(RCS.dg),'r-','LineWidth',2);
% hold on
% plot((180/pi)*(theta_obser+pi),10*log10(RCS.pod),'b-','LineWidth',2)
% legend('DGTD','POD-DGTD')
% xlabel('\theta ');
% ylabel('RCS/dBsm');
% grid on

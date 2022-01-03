function PODsolution = poddgsolution_dft(ANNsolution)
%% DFT
global parameter;
%% using fft
T = parameter.dt;
Fs = 1/T;
N = size(ANNsolution,2);
NFFT = 2^nextpow2(N);
% Fx = 0:Fs/(N-1):Fs;
Fx = (0:N-1)*(Fs/N);
Fn = parameter.freq/parameter.c0;
fftf = fft(ANNsolution');
% PODsolution = fftf(2,:)'/(N/2);
PODsolution(:,1) = interp1(Fx',fftf,Fn)/(N/2);
%%
% x = GPRsolution;
% M = size(x,1);
% N = size(x,2);
% omegN = exp(-2*pi*(1i)/N);
% PODsolution = zeros(M,N);
% omegNk = 
% for k = 1:N
%     omegNk = omegN.^(((1:N)-1)*(k - 1));
% end
% 
% for ii = 1:M
%     for k = 1:N
%         PODsolution(ii,k) = sum(x(ii,1:N).*);
%     end
% end
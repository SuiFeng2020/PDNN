%% 计算误差
Np = 81;
Nt = 263;
% 计算数值解Uh 与 投影解Upro 的误差
Errpro_Ez = 
for i = 1:Np
    for j = 1:Nt
        Phi_Uh = timeparameterPOD.Basis.Ez' * Snapshots(i).Ez(:,j);
        Upro = timeparameterPOD.Basis.Ez * Phi_Uh;
        Err = Snapshots(i).Ez(:,j) - Upro;
    end
end

        
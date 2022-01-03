function traintime = gettraintime(parameter)
time = 0;
deltat = parameter.dt;
tt = 1;
while time < parameter.tmax*parameter.c0
    % Increment time
    time  = time + deltat;
    if time >= (parameter.tmax*parameter.freq - 1)
        traintime(tt) = time;
        tt = tt + 1;
    end
end
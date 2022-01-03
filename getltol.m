function tol = getltol(l)
% SVD truncation tolerance
if l <= 5 
        tol = 0.01;
    elseif l >5 && l <= 10
        tol = 0.05;
    elseif l >10 && l <= 20
        tol = 0.1;
    elseif l >20 && l <= 30
        tol = 0.2;
    elseif l >30 && l <= 40
        tol = 0.3;
    elseif l >40 && l <= 55
        tol = 0.4;
    else
        tol = 0.5;
end
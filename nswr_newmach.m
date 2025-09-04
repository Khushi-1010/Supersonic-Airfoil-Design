function M = nswr_newmach(M1)
    gamma = 1.4;  % For air
    M = sqrt((1+(gamma-1)*M1^2/2)/((gamma*M1^2)-(gamma-1)/2));
end
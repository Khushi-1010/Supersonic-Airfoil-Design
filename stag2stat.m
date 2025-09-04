function P0P = stag2stat(M)
    gamma = 1.4;
    P0P = (1+(gamma-1)*M^2/2)^(gamma/(gamma-1));
end
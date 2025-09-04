function p2p1 = nswr_pressure(M)
    gamma = 1.4;
    p2p1 = 1 + ((2*gamma/(gamma+1))*(M^2-1));
end
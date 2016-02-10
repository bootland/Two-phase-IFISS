function [tpve,tpde] = bottomlefttriangle(nx,visc1,visc2,dens1,dens2,qmethod)

if qmethod == 0
    tpn = (nx-1);
elseif qmethod == 2
    tpn = (nx-1)/2;
end
TPE = zeros(tpn,tpn);
for k = 1:tpn-1
    TPE(tpn-k,1:k) = 1;
end
tpve = visc1 + (visc2-visc1)*TPE(:);
tpde = dens1 + (dens2-dens1)*TPE(:);

end
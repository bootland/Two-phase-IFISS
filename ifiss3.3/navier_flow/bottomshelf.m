function [tpve,tpde] = bottomshelf(nx,visc1,visc2,dens1,dens2,qmethod)

if qmethod == 0
    tpn = (nx-1);
elseif qmethod == 2
    tpn = (nx-1)/2;
end
TPE = zeros(tpn,tpn);
TPE(:,1:tpn/2) = 1;
tpve = visc1 + (visc2-visc1)*TPE(:);
tpde = dens1 + (dens2-dens1)*TPE(:);

end
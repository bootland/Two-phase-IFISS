function [tpv,tpd,tpve,tpvh] = toprightblock(nv,visc1,visc2,dens1,dens2,qmethod)

tpn = sqrt(nv)-1;
TP = zeros(tpn+1,tpn+1);
TP(tpn/2+1:end,tpn/2+1:end) = 1;
tpv = visc1 + (visc2-visc1)*TP;
if qmethod < 2
    tpd = dens1 + (dens2-dens1)*TP(:);
    tpvh = visc1 + (visc2-visc1)*TP(:);
else
    TPH = TP(1:2:end,1:2:end);
    tpd = dens1 + (dens2-dens1)*TPH(:);
    tpvh = visc1 + (visc2-visc1)*TPH(:);
end

tpve = zeros(tpn,tpn);
for j = 1:tpn
    for k = 1:tpn
        tpve(j,k) = geomean([tpv(j,k:k+1),tpv(j+1,k:k+1)]);
    end
end
tpv = tpv(:);
tpve = tpve(:);

end
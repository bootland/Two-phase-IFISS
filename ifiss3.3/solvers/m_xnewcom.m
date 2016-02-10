function y = m_xnewcom(x_it,aparams,mparams)
%M_XBFBT ideal least squares commutator preconditioner
%   y = m_xbfbt(x_it,aparams,mparams);
%   input
%          x_it         operand for preconditioning operator
%          aparams      structure defining coefficient matrix
%          mparams      structure defining preconditioning matrix
%   output
%          y            result of preconditioning operation
%
%   IFISS function: HCE; 15 March 2005.
% Copyright (c) 2005 D.J. Silvester, H.C. Elman, A. Ramage

nv = length(aparams.F);
%nu = nv/2;
np = size(aparams.B,1);
% Gdiag=spdiags(diag(mparams.G),0,nv,nv);
Mpdiag=spdiags(diag(mparams.Mp),0,np,np);

rv=x_it(1:nv); rp=x_it(nv+1:nv+np);

%% pressure solve
% xB = (Gdiag\aparams.B')';
xB = Mpdiag\aparams.B;
% BBt = aparams.B*xB';
   
if mparams.domain==1,
%    n_null = mparams.n_null;
%    minor = [1:n_null-1,n_null+1:np]';
%    rp1 = zeros(np,1);
%    rp1 = (aparams.F)*((mparams.A)\(xB*rp));
%    rp2 =  (mparams.A)\rp1;
%    zp = zeros(np,1);
%    zp = -(xB')*rp2;
    zp = - (xB)*(aparams.F*((mparams.A)\(xB'*rp)));
%    rp1 = zeros(np,1);
%    rp1(minor) = Mpdiag(minor,minor)\rp(minor);
%    rp2 = aparams.B*((mparams.A)\(aparams.F*((mparams.A)\(aparams.B'*rp1))));
%    zp = zeros(np,1);
%    zp(minor) = - Mpdiag(minor,minor)\rp2(minor);
else    
%    zp = -BBt \ (xB*(aparams.F*(xB'*(BBt\rp))));
   zp = - (xB')*(aparams.F*((mparams.A)\(xB*rp)));
end 

%% velocity solve
rv = rv-(aparams.B')*zp;
zv = aparams.F \ rv;
y = [zv;zp];

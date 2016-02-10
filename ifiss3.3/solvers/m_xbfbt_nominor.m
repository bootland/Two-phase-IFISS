function y = m_xbfbt_nominor(x_it,aparams,mparams)
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
Gdiag=spdiags(diag(mparams.G),0,nv,nv);

rv=x_it(1:nv); rp=x_it(nv+1:nv+np);

%% pressure solve
xB = (Gdiag\aparams.B')';
BBt = aparams.B*xB';
   
zp = -BBt \ (xB*(aparams.F*(xB'*(BBt\rp))));

%% velocity solve
rv = rv-(aparams.B')*zp;
zv = aparams.F \ rv;
y = [zv;zp];

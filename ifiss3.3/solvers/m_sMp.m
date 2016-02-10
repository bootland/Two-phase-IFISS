function y = m_sMp(x_it,aparams,mparams)
%M_SMP ideal scaled mass matrix preconditioner
%   y = m_sMp(x_it,aparams,mparams)
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
np = size(aparams.B,1);

rv=x_it(1:nv); rp=x_it(nv+1:nv+np);

sMp = (mparams.Mp)./(mparams.viscosity);

%% pressure solve
if mparams.domain==1,
   n_null = mparams.n_null;
   minor = [1:n_null-1,n_null+1:np]';
   yp = zeros(np,1);
   yp(minor) = sMp(minor,minor)\rp(minor);
%    yp = sMp\rp;
   zp = -yp;
else
   zp = -sMp\rp;
end

%% velocity solve
rv = rv-(aparams.B')*zp;
zv = aparams.F \ rv;
y = [zv;zp];
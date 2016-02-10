function y = m_Mp(x_it,aparams,mparams)
%M_SMP2P ideal scaled mass matrix preconditioner
%   y = m_sMp2p(x_it,aparams,mparams)
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

%% pressure solve
if mparams.domain==1,
   n_null = mparams.n_null;
   minor = [1:n_null-1,n_null+1:np]';
   zp = zeros(np,1);
%    zp(minor) = mparams.Mp(minor,minor)\rp(minor);
   zp(minor) = - mparams.Mp(minor,minor)\rp(minor);
else
%    zp = mparams.Mp\rp;
   zp = -mparams.Mp\rp;
end

%% velocity solve
rv = rv-(aparams.B')*zp;
zv = aparams.F \ rv;
y = [zv;zp];
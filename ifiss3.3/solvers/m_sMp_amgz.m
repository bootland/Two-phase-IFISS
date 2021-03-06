function y = m_sMp_amgz(x_it,aparams,mparams)
%M_SMP_AMGZ AMG iterated scaled mass matrix preconditioner
%   y = m_sMp_amgz(x_it,aparams,mparams)
%   input
%          x_it         operand for preconditioning operator
%          aparams      structure defining coefficient matrix
%          mparams      structure defining preconditioning matrix
%   output
%          y            result of preconditioning operation
%
%   IFISS function: HCE; 15 March 2005.
% Copyright (c) 2005 D.J. Silvester, H.C. Elman, A. Ramage

global  amg_gridA amg_smootherA amg_gridF amg_smootherF

nv = length(aparams.F);
np = size(aparams.B,1);

rv=x_it(1:nv); rp=x_it(nv+1:nv+np);

%% pressure solve
if mparams.domain==1,
   n_null = mparams.n_null;
   minor = [1:n_null-1,n_null+1:np]';
   zp = zeros(np,1);
   zp(minor) = - amg_v_cycle(rp(minor), amg_gridA,  amg_smootherA);
else
   zp = -amg_v_cycle(rp, amg_gridA,  amg_smootherA);
end

%% velocity solve
rv = rv-(aparams.B')*zp;
zv = amg_v_cycle(rv, amg_gridF,  amg_smootherF);
y = [zv;zp];
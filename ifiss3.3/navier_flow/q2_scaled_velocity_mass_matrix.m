function [r] = q2_scaled_velocity_mass_matrix(xy,~,mv,~,tpse)
%Q2_SCALED_VELOCITY_MASS_MATRIX Q2 scaled velocity mass matrix generator
%   [r] = q2_scaled_velocity_mass_matrix(xy,xyp,mv,mp,tpse);
%   input
%          xy         Q2 nodal coordinate vector
%          xyp        Q1 nodal coordinate vector
%          mv         Q2 element mapping matrix
%          mp         Q1 element mapping matrix
%          tpse       scaling vector, two-phase elementwise scaling
%   output
%          r          scaled Q2 vector mass matrix
%
%   Natural boundary conditions apply. Dirichlet conditions
%   must be explicitly enforced by calling function flowbc.
%   IFISS function: DJS; 6 March 2005.
% Copyright (c) 2005 D.J. Silvester, H.C. Elman, A. Ramage
nngpt=9;
x=xy(:,1); y=xy(:,2);
% xp=xyp(:,1); yp=xyp(:,2);
nvtx=length(x); nu=2*nvtx; % np=length(xp);
nel=length(mv(:,1));
fprintf('setting up Q2 scaled velocity mass matrix... \n')
%
% initialise global matrices
r = sparse(nu,nu);
%
% Gauss point integration rules
if (nngpt==4)        % 2x2 Gauss points
    gpt=1.0e0/sqrt(3.0e0);
    s(1) = -gpt; t(1) = -gpt; wt(1)=1;
    s(2) =  gpt; t(2) = -gpt; wt(2)=1;
    s(3) =  gpt; t(3) =  gpt; wt(3)=1;
    s(4) = -gpt; t(4) =  gpt; wt(4)=1;
elseif (nngpt==9)   % 3x3 Gauss points
    gpt=sqrt(0.6);
    s(1) = -gpt; t(1) = -gpt; wt(1)=25/81;
    s(2) =  gpt; t(2) = -gpt; wt(2)=25/81;
    s(3) =  gpt; t(3) =  gpt; wt(3)=25/81;
    s(4) = -gpt; t(4) =  gpt; wt(4)=25/81;
    s(5) =  0.0; t(5) = -gpt; wt(5)=40/81;
    s(6) =  gpt; t(6) =  0.0; wt(6)=40/81;
    s(7) =  0.0; t(7) =  gpt; wt(7)=40/81;
    s(8) = -gpt; t(8) =  0.0; wt(8)=40/81;
    s(9) =  0.0; t(9) =  0.0; wt(9)=64/81;
else
    error('Check Gauss point integration specification')
end
%
% inner loop over elements
xl_v = zeros(nel,4);
yl_v = zeros(nel,4);
for ivtx = 1:4
    xl_v(:,ivtx) = x(mv(:,ivtx));
    yl_v(:,ivtx) = y(mv(:,ivtx));
end
re = zeros(nel,9,9);
%
% loop over Gauss points
for igpt = 1:nngpt
    sigpt=s(igpt);
    tigpt=t(igpt);
    wght=wt(igpt);
    %  evaluate derivatives etc
    [jac,~,~,~,~] = deriv(sigpt,tigpt,xl_v,yl_v);
    [psi,~,~] = qderiv(sigpt,tigpt,xl_v,yl_v);
    for j = 1:9
        for i = 1:9
            re(:,i,j)  = re(:,i,j)  + wght*tpse.*psi(:,i).*psi(:,j).*jac(:);
        end
    end
    % end of Gauss point loop
end
%
% element assembly into global matrices
for krow=1:9
    nrow=mv(:,krow);
    for kcol=1:9
        ncol=mv(:,kcol);
        r = r + sparse(nrow,ncol,re(:,krow,kcol),nu,nu);
        r = r + sparse(nrow+nvtx,ncol+nvtx,re(:,krow,kcol),nu,nu);
    end
end
%
fprintf('done\n')
return

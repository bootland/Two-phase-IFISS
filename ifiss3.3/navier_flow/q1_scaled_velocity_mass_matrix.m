function [r] = q1_scaled_velocity_mass_matrix(xy,ev,tpse)
%Q1_SCALED_VELOCITY_MASS_MATRIX Q1 scaled velocity mass matrix generator
%   [r] = q1_scaled_velocity_mass_matrix(xy,ev,tpse);
%   input
%          xy         Q2 nodal coordinate vector 
%          ev         element mapping matrix
%          tpse       scaling vector, two-phase elementwise scaling
%   output
%          r          scaled Q1 vector mass matrix
%
%   Natural boundary conditions apply. Dirichlet conditions
%   must be explicitly enforced by calling function flowbc.
%   IFISS function: DJS; 8 March 2005.
% Copyright (c) 2005 D.J. Silvester, H.C. Elman, A. Ramage 
nngpt=4; 
x=xy(:,1); y=xy(:,2);
% xp=xy(:,1); yp=xy(:,2);
nvtx=length(x); nu=2*nvtx; % np=length(xp); 
nel=length(ev(:,1));
fprintf('setting up Q1 scaled velocity mass matrix... \n')
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
elseif (nngpt==1)   % 1x1 Gauss point
   s(1) =    0; t(1) =    0; wt(1)=4;
else
   error('Check Gauss point integration specification')
end
%
% inner loop over elements
xl_v = zeros(nel,4);
yl_v = zeros(nel,4);
for ivtx = 1:4
   xl_v(:,ivtx) = x(ev(:,ivtx));
   yl_v(:,ivtx) = y(ev(:,ivtx)); 
end
re = zeros(nel,9,9);
% 
% loop over Gauss points
for igpt = 1:nngpt
   sigpt=s(igpt);
   tigpt=t(igpt);
   wght=wt(igpt);
%  evaluate derivatives etc
   [jac,~,phi,~,~] = deriv(sigpt,tigpt,xl_v,yl_v);
   for j = 1:4
      for i = 1:4
         re(:,i,j)  = re(:,i,j)  + wght*tpse.*phi(:,i).*phi(:,j).*jac(:);
      end
   end
% end of Gauss point loop
end  
%
% element assembly into global matrices
for krow=1:4
   nrow=ev(:,krow);	 
   for kcol=1:4
      ncol=ev(:,kcol);	  
      r = r + sparse(nrow,ncol,re(:,krow,kcol),nu,nu);
      r = r + sparse(nrow+nvtx,ncol+nvtx,re(:,krow,kcol),nu,nu);
   end
end
%
fprintf('done\n')
return

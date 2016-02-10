function [m] = q1_scaled_pressure_mass_matrix(xy,xyp,mv,mp,tpse)
%Q1_SCALED_PRESSURE_MASS_MATRIX Q1 scaled pressure mass matrix generator
%   [m] = q1_scaled_pressure_mass_matrix(xy,xyp,mv,mp,tpse);
%   input
%          xy         Q2 nodal coordinate vector
%          xyp        Q1 element coordinate vector
%          mv|ev      Q2|Q1 element mapping matrix
%          mp|ev      Q1|Q1 element mapping matrix
%          tpse       scaling vector, two-phase elementwise scaling
%   output
%          m          scaled Q1 mass matrix 
%
%   Natural boundary conditions apply. Dirichlet conditions
%   must be explicitly enforced by calling function flowbc.
%   IFISS function: DJS; 8 March 2005.
% Copyright (c) 2005 D.J. Silvester, H.C. Elman, A. Ramage 
nngpt=4; 
x=xy(:,1); y=xy(:,2);
xp=xyp(:,1); % yp=xyp(:,2);
np=length(xp); % nvtx=length(x); nu=2*nvtx;
nel=length(mv(:,1));
fprintf('setting up Q1 scaled pressure mass matrix... \n')
%
% initialise global matrices
m = sparse(np,np);
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
   xl_v(:,ivtx) = x(mv(:,ivtx));
   yl_v(:,ivtx) = y(mv(:,ivtx)); 
end
mpe = zeros(nel,4,4);
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
         mpe(:,i,j) = mpe(:,i,j) + wght*tpse.*phi(:,i) .*phi(:,j) .*jac(:);
      end
   end
% end of Gauss point loop
end  
%
% element assembly into global matrices
for krow=1:4
   nrow=mp(:,krow);	 
   for kcol=1:4
      ncol=mp(:,kcol);	  
      m = m + sparse(nrow,ncol,mpe(:,krow,kcol),np,np);
   end
end
%
fprintf('done\n')
return

function [a,b,bbx,bby,muc,rhor,f,g,l] = stokes_q1q1_two_phase_on_elements_plotting(xy,ev,tpve,tpde)
%STOKES_Q1Q1_TWO_PHASE_ON_ELEMENTS_PLOTTING Q1-Q1 matrix generator
%   [A,B,Bx,By,muC,rhoG,f,g,L] = stokes_q1q1_two_phase_on_elements_plotting(xy,ev,tpve,tpde);
%   input
%          xy         Q2 nodal coordinate vector 
%          ev         element mapping matrix
%          tpve       vector of the two-phase viscosity on elements
%          tpde       vector of the two-phase density on elements
%   output
%          A          viscosity scaled Q1 vector diffusion matrix
%          B          Q1-Q1 divergence matrix
%          Bx         Q1 x-derivative matrix
%          By         Q1 y-derivative matrix
%          muC        viscosity scaled pressure stabilization matrix
%          rhoG       density scaled Q2 vector mass matrix
%          f          velocity rhs vector
%          g          pressure rhs vector
%          L          unscaled Q2 vector diffusion matrix (for plotting)
%
%   Natural boundary conditions apply. Dirichlet conditions
%   must be explicitly enforced by calling function flowbc.
%   IFISS function: DJS; 8 March 2005.
% Copyright (c) 2005 D.J. Silvester, H.C. Elman, A. Ramage 
nngpt=4; 
x=xy(:,1); y=xy(:,2);
xp=xy(:,1); % yp=xy(:,2);
nvtx=length(x); nu=2*nvtx; np=length(xp); 
nel=length(ev(:,1)); % mp=[1:nel]';
% lx=max(x)-min(x); ly=max(y)-min(y);
% hx=max(diff(x)); hy=max(diff(y));
fprintf('setting up Q1-Q1 matrices...  ')
%
% initialise global matrices
a = sparse(nu,nu);
l = sparse(nu,nu);
rhor = sparse(nu,nu);
bbx = sparse(nvtx,nvtx);
bby = sparse(nvtx,nvtx);
bx = sparse(np,nvtx);
by = sparse(np,nvtx);
% b = sparse(np,nu);
muc = sparse(np,np);
f = zeros(nu,1);
g = zeros(np,1);
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
ae = zeros(nel,4,4);
le = zeros(nel,4,4);
rhore = zeros(nel,9,9);
bbxe = zeros(nel,4,4);
bbye = zeros(nel,4,4);
bxe = zeros(nel,4,4);
bye = zeros(nel,4,4);
% ge = zeros(nel,4);
% 
% loop over Gauss points
for igpt = 1:nngpt
   sigpt=s(igpt);
   tigpt=t(igpt);
   wght=wt(igpt);
%  evaluate derivatives etc
   [jac,invjac,phi,dphidx,dphidy] = deriv(sigpt,tigpt,xl_v,yl_v);
   for j = 1:4
      for i = 1:4
         ae(:,i,j)  = ae(:,i,j)  + wght*tpve.*dphidx(:,i).*dphidx(:,j).*invjac(:);
         ae(:,i,j)  = ae(:,i,j)  + wght*tpve.*dphidy(:,i).*dphidy(:,j).*invjac(:);
         le(:,i,j)  = le(:,i,j)  + wght*dphidx(:,i).*dphidx(:,j).*invjac(:);
         le(:,i,j)  = le(:,i,j)  + wght*dphidy(:,i).*dphidy(:,j).*invjac(:);
         rhore(:,i,j)  = rhore(:,i,j)  + wght*tpde.*phi(:,i).*phi(:,j).*jac(:);
         bbxe(:,i,j) = bbxe(:,i,j) - wght*phi(:,i) .*dphidx(:,j);             
         bbye(:,i,j) = bbye(:,i,j) - wght*phi(:,i) .*dphidy(:,j);   
         bxe(:,i,j) = bxe(:,i,j) - wght* dphidx(:,j).*phi(:,i);
         bye(:,i,j) = bye(:,i,j) - wght* dphidy(:,j).*phi(:,i);
      end
   end
% end of Gauss point loop
end  
%
% element assembly into global matrices
% component velocity matrices ...    
for krow=1:4
   nrow=ev(:,krow);	 
   for kcol=1:4
      ncol=ev(:,kcol);	  
      a = a + sparse(nrow,ncol,ae(:,krow,kcol),nu,nu);
      a = a + sparse(nrow+nvtx,ncol+nvtx,ae(:,krow,kcol),nu,nu);
      l = l + sparse(nrow,ncol,le(:,krow,kcol),nu,nu);
      l = l + sparse(nrow+nvtx,ncol+nvtx,le(:,krow,kcol),nu,nu);
      rhor = rhor + sparse(nrow,ncol,rhore(:,krow,kcol),nu,nu);
      rhor = rhor + sparse(nrow+nvtx,ncol+nvtx,rhore(:,krow,kcol),nu,nu);
      bbx = bbx + sparse(nrow,ncol,bbxe(:,krow,kcol),nvtx,nvtx);
      bby = bby + sparse(nrow,ncol,bbye(:,krow,kcol),nvtx,nvtx);
      bx = bx + sparse(ncol,nrow,bxe(:,kcol,krow),nvtx,nvtx);
      by = by + sparse(ncol,nrow,bye(:,kcol,krow),nvtx,nvtx);
   end
end
%
% vector velocity matrices ...
b = [bx,by];
%   
%%
% set up Bochev pressure stabilisation matrix
mucpe = zeros(nel,4,4);
% loop over Gauss points
for igpt = 1:nngpt
   sigpt=s(igpt);
   tigpt=t(igpt);
   wght=wt(igpt);
   [jac,~,phi,~,~] = deriv(sigpt,tigpt,xl_v,yl_v);
   for j = 1:4
      for i = 1:4
%%%% NW improved version
         mucpe(:,i,j) = mucpe(:,i,j) + (wght./tpve).*(phi(:,i)-0.25) .*(phi(:,j)-0.25) .*jac(:);
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
      muc = muc + sparse(nrow,ncol,mucpe(:,krow,kcol),nvtx,nvtx);
   end
end
fprintf('done\n')
return

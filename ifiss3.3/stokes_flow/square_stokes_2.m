%SQUARE_STOKES set up flow problem in unit square domain 
%   IFISS scriptfile: DJS; 27 May 2012.
% Copyright (c) 2005 D.J. Silvester, H.C. Elman, A. Ramage 
clear variables
%% define geometry
pde=3; domain=1; enclosed=1;
cavity_domain
load cavity_grid.mat
%
%% set up matrices
load('testrun.mat')
if testrun ~= 1
    q_in=default('Q1-Q1/Q1-P0/Q2-Q1/Q2-P1: 1/2/3/4? (default Q2-Q1)',3);
else
    q_in = 3;
end
qmethod=q_in-1;
if (qmethod==1 || qmethod==3)
    error('Element space choice not yet supported for two-phase problems.')
end
if testrun ~= 1
    phstr=default('Phase structure: middle square/bottom shelf/bottom left triangle 1/2/3',1);
    visc1=default('kinematic viscosity parameter of "outer" phase 1 (default 2/1000)',2/1000);
    visc2=default('kinematic viscosity parameter of "inner" phase 2 (default 2/10)',2/10);
    dens1=default('density parameter of "outer" phase 1 (default 1)',1);
    dens2=default('density parameter of "inner" phase 2 (default 100)',100);
    nx = length(x);
    if phstr==1
        [tpkve,tpde] = middlehalfblock(nx,visc1,visc2,dens1,dens2,qmethod);
    elseif phstr==2
        [tpkve,tpde] = bottomshelf(nx,visc1,visc2,dens1,dens2,qmethod);
    elseif phstr==3
        [tpkve,tpde] = bottomlefttriangle(nx,visc1,visc2,dens1,dens2,qmethod);
    else
        error('Not a valid choice of phase structure')
    end
else
    phstr = 1;
    visc1 = 2/Re;
    visc2 = visc1*viscd;
    dens1 = 1;
    dens2 = densd;
    if phstr==1
        [tpkve,tpde] = middlehalfblock(length(x),visc1,visc2,dens1,dens2,qmethod);
    elseif phstr==2
        [tpkve,tpde] = bottomshelf(length(x),visc1,visc2,dens1,dens2,qmethod);
    elseif phstr==3
        [tpkve,tpde] = bottomlefttriangle(length(x),visc1,visc2,dens1,dens2,qmethod);
    end
end
tpve = tpkve.*tpde;
if qmethod==2,
   [x,y,xy,xyp,mp,map] = q2q1gridx(x,y,xy,mv,bound);
   if testrun ~= 1
       [A,B,Bx,By,rhoG,f,g,L] = stokes_q2q1_two_phase_on_elements_plotting(xy,xyp,mv,mp,tpve,tpde);
   else
       [A,B,Bx,By,rhoG,f,g] = stokes_q2q1_two_phase_on_elements(xy,xyp,mv,mp,tpve,tpde);
   end
elseif qmethod==3,
   [x,y,xy,xyp,ee] = q2p1grid(x,y,xy,mv,bound);
   [A,B,Q,G,Bx,By,f,g] = stokes_q2p1(xy,xyp,mv);
elseif qmethod==8,
   [x,y,xy,xyp,ee] = q2p1grid(x,y,xy,mv,bound);
   [A,BB,QQ,G,Bx,By,f,gg] = stokes_q2p1(xy,xyp,mv);
   nnp=length(gg); ppk=[1:3:nnp];
   B=BB(ppk,:); Q=QQ(ppk,ppk); g=gg(ppk);
   fprintf('Reduction to Q2-P0 Stokes system done.\n')
elseif qmethod==0, 
   [ev,ee,ebound,xyp] = q1q1grid(x,y,xy,mv,bound,mbound);
   [~,elinds] = sort(ev(:,1)); [~,ninds] = sort(elinds); tpve = tpve(ninds); tpde = tpde(ninds);
   if testrun ~= 1
       [A,B,Bx,By,muC,rhoG,f,g,L] = stokes_q1q1_two_phase_on_elements_plotting(xy,ev,tpve,tpde);
   else
       [A,B,Bx,By,muC,rhoG,f,g] = stokes_q1q1_two_phase_on_elements(xy,ev,tpve,tpde);
   end
elseif qmethod==1 
   [ev,ee,ebound,xyp] = q1p0grid(x,y,xy,mv,bound,mbound);
   [A,B,Q,C,G,Bx,By,f,g] = stokes_q1p0(xy,xyp,mv,ev);
end
gohome
cd datafiles
save square_stokes_nobc.mat pde domain qmethod grid_type A B f g xy xyp mbound bound x y 
save square_stokes_nobc.mat Bx By bndxy bnde obs -append
save square_stokes_nobc.mat tpkve tpve tpde visc1 visc2 dens1 dens2 -append
if qmethod==1 
   save square_stokes_nobc.mat C G ev ee ebound -append
elseif qmethod==0 
   save square_stokes_nobc.mat muC rhoG ev ee ebound mv enclosed -append
elseif qmethod==2
   save square_stokes_nobc.mat mv mp rhoG map -append
else
   save square_stokes_nobc.mat mv ee G -append
end
if (qmethod==0 || qmethod==2) && testrun~=1
   save square_stokes_nobc.mat L -append
end
fprintf('system matrices saved in square_stokes_nobc.mat ...\n')

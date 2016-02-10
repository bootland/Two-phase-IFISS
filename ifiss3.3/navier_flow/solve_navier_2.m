%SOLVE_NAVIER solve singular Navier-Stokes problem
%   IFISS scriptfile: DJS; 29 September 2013.
% Copyright (c) 2005 D.J. Silvester, H.C. Elman, A. Ramage 
clear variables
load('testrun.mat')
fprintf('Enclosed flow problem ...\n')

if testrun ~= 1
    deltat=default('Time step (choose 0 for steady problem, this is default)',0);
    nlmethod=default('Picard/Newton/hybrid linearization 1/2/3 (default Picard)',1);
    nlmethod=nlmethod-1;
    if nlmethod==0,
        maxit_p=default('number of Picard iterations (default 1)',1);
        maxit_n=0;
    elseif nlmethod==1,
        maxit_p=0;
        maxit_n=default('number of Newton iterations (default 6)',6);
    else
        maxit_p=default('number of Picard iterations (default 2)',2);
        maxit_n=default('number of Newton iterations (default 4)',4);
    end
    tol_nl=default('nonlinear tolerance (default 1.e-8)',1.e-8);
else
    deltat = timestep;
    nlmethod = 0;
    maxit_p = 1;
    maxit_n = 0;
    tol_nl = 1.e-8;
end
%
%
%% initialize for nonlinear iteration: compute Stokes solution
%% load assembled matrices
gohome
cd datafiles
load square_stokes_nobc.mat
%
fprintf('stokes system ...\n')

%% Add in mass matrix term for nonsteady probelms with timestep deltat

if deltat ~= 0
    A = A + rhoG/deltat;
end

%% boundary conditions
[Ast,Bst,fst,gst] = flowbc(A,B,f,g,xy,bound);
nlres0_norm = norm([fst;gst]);
%
nv=length(fst)/2; np=length(gst);

if qmethod>1,
   beta=0;
% xst=[Ast,Bst';Bst,sparse(np,np)]\[fst;gst];
%----------------------------------- stabilized version
xstz=[Ast,Bst',zeros(2*nv,1);Bst,sparse(np,np),ones(np,1)/np; ...
       zeros(1,2*nv),ones(1,np)/np,zeros(1,1)]\[fst;gst;0];
xst=xstz(1:end-1); multiplier=xstz(end);
elseif qmethod==1
	  beta=1/4;     % default parameter
%  xst=[Ast,Bst';Bst,-beta*C]\[fst;gst];
%----------------------------------- stabilized version
xstz=[Ast,Bst',zeros(2*nv,1);Bst,-beta*C,ones(np,1)/np; ...
      zeros(1,2*nv),ones(1,np)/np,zeros(1,1)]\[fst;gst;0];
xst=xstz(1:end-1); multiplier=xstz(end);
elseif qmethod==0
   fprintf('computing pressure stabilized solution...\n')
   beta=1;
%  xst=[Ast,Bst';Bst,-beta*C]\[fst;gst];
%----------------------------------- stabilized version
% % xstz=[Ast,Bst',zeros(2*nv,1);Bst,-beta*C,ones(np,1)/np; ...
% %       zeros(1,2*nv),ones(1,np)/np,zeros(1,1)]\[fst;gst;0];
xstz=[Ast,Bst',zeros(2*nv,1);Bst,-beta*muC,ones(np,1)/np; ...
      zeros(1,2*nv),ones(1,np)/np,zeros(1,1)]\[fst;gst;0];
xst=xstz(1:end-1); multiplier=xstz(end);
end

% compute residual of Stokes solution
if qmethod>1
   N = navier_q2(xy,mv,tpde,xst);
elseif qmethod<=1,
%    nubeta=beta/viscosity;
%    nubeta=tpdi*spdiags(beta./tpv,0,tplen,tplen);
% %    nubeta = beta/sqrt(visc1*visc2);
   N = navier_q1(xy,ev,tpde,xst);
end
% Anst = viscosity*A + [N, sparse(nv,nv); sparse(nv,nv), N];
% Anst = tpvd*A + [N, sparse(nv,nv); sparse(nv,nv), N];
Anst = A + [N, sparse(nv,nv); sparse(nv,nv), N];
[Anst,Bst,fst,gst] = flowbc(Anst,B,f,g,xy,bound);
% Bst = tpddi*Bst; gst = tpddi*gst;
if     qmethod>1, nlres = [Anst,Bst';Bst,sparse(np,np)]*xst-[fst;gst];
elseif qmethod<=1, % % nlres = [Anst,Bst';Bst,-nubeta*C]*xst-[fst;gst];
    nlres = [Anst,Bst';Bst,-beta*muC]*xst-[fst;gst];
end
nlres_norm  = norm(nlres);
%%% plot solution
if testrun ~= 1
    spc=default('uniform/nonuniform streamlines 1/2 (default nonuniform)',2);
else
    spc = 2;
end
%contourn = default('number of contour lines (default 50)',50);
if testrun ~= 1
%     flowplot(qmethod,xst,By,Bx,A,xy,xyp,x,y,bound,spc,33);
    flowplot(qmethod,xst,By,Bx,L,xy,xyp,x,y,bound,spc,33);
    %flowplot09(qmethod,xst,By,Bx,A,xy,xyp,x,y,bound,bndxy,bnde,obs,contourn,spc,133)
    pause(1)
end
fprintf('\n\ninitial nonlinear residual is %e ',nlres0_norm)
fprintf('\nStokes solution residual is %e\n', nlres_norm)
flowsol = xst;
%
%
pde=4;
it_p = 0;
%
% nonlinear iteration 
%% Picard startup step
while nlres_norm>nlres0_norm*tol_nl && it_p<maxit_p,
  %nlres = nlres - [zeros(2*nv,1);(sum(nlres(2*nv+1:2*nv+np))/np)*ones(np,1)];
   it_p = it_p+1;
   fprintf('\nPicard iteration number %g \n',it_p),
% compute Picard correction and update solution
   if     qmethod>1,
%  dxns = -[Anst,Bst';Bst,sparse(np,np)]\nlres;
%----------------------------------- stabilized version
   dxnsz= -[Anst,Bst',zeros(2*nv,1);Bst,sparse(np,np),ones(np,1)/np; ...
                zeros(1,2*nv),ones(1,np)/np,zeros(1,1)]\[nlres;0];
   dxns=dxnsz(1:end-1); multiplier=dxnsz(end);
   elseif qmethod<=1,
%  dxns = -[Anst,Bst';Bst,-nubeta*C]\nlres;
%----------------------------------- stabilized version
% %    dxnsz= -[Anst,Bst',zeros(2*nv,1);Bst,-nubeta*C,ones(np,1)/np; ...
% %             zeros(1,2*nv),ones(1,np)/np,zeros(1,1)]\[nlres;0];
   dxnsz= -[Anst,Bst',zeros(2*nv,1);Bst,-beta*muC,ones(np,1)/np; ...
            zeros(1,2*nv),ones(1,np)/np,zeros(1,1)]\[nlres;0];
   dxns=dxnsz(1:end-1); multiplier=dxnsz(end);
   end
   xns = flowsol + dxns;
% compute residual of new solution
   if     qmethod>1,  N = navier_q2(xy,mv,tpde,xns);
   elseif qmethod<=1, N = navier_q1(xy,ev,tpde,xns);
   end
%    Anst = tpvd*A + [N, sparse(nv,nv); sparse(nv,nv), N];
   Anst = A + [N, sparse(nv,nv); sparse(nv,nv), N];
   [Anst,Bst,fst,gst] = flowbc(Anst,B,f,g,xy,bound);
%    Bst = tpddi*Bst; gst = tpddi*gst;
   if     qmethod>1,  nlres = [Anst,Bst';Bst,sparse(np,np)]*xns-[fst;gst];
   elseif qmethod<=1, % % nlres = [Anst,Bst';Bst,-nubeta*C]*xns-[fst;gst];
       nlres = [Anst,Bst';Bst,-beta*muC]*xns-[fst;gst];
   end
   nlres_norm = norm(nlres);
   nnv=length(fst); soldiff=norm(xns(1:nnv)-flowsol(1:nnv));
   fprintf('nonlinear residual is %e',nlres_norm)
   fprintf('\n   velocity change is %e\n',soldiff)
   if testrun ~= 1
       % plot solution
%        flowplot(qmethod,xns,By,Bx,A,xy,xyp,x,y,bound,spc,66); drawnow;
       flowplot(qmethod,xns,By,Bx,L,xy,xyp,x,y,bound,spc,66); drawnow;
       %  flowplot09(qmethod,xns,By,Bx,A,xy,xyp,x,y,bound,bndxy,bnde,obs,contourn,spc,166);drawnow
       pause(1)
   end
   flowsol = xns;
%% end of Picard iteration loop
end
%%
%
it_nl = it_p;
it_n = 0;
%% Newton iteration loop
while (nlres_norm > nlres0_norm*tol_nl) && (it_nl < maxit_p + maxit_n),
  %nlres = nlres - [zeros(2*nv,1);(sum(nlres(2*nv+1:2*nv+np))/np)*ones(np,1)];
   it_n = it_n+1;
   it_nl = it_nl+1;
   fprintf('\nNewton iteration number %g \n',it_n),
% compute Jacobian of current solution
   if     qmethod>1,  [Nxx,Nxy,Nyx,Nyy] = newton_q2(xy,mv,flowsol);
   elseif qmethod<=1, [Nxx,Nxy,Nyx,Nyy] = newton_q1(xy,ev,flowsol);
   end
%    J = tpvd*A + [N + Nxx, Nxy; Nyx, N + Nyy];
   J = A + [N + Nxx, Nxy; Nyx, N + Nyy];
   Jnst = newtonbc(J,xy,bound); 
% compute Newton correction and update solution
   if qmethod>1,
%  dxns = -[Jnst,Bst';Bst,sparse(np,np)]\nlres;
%----------------------------------- stabilized version
   dxnsz= -[Jnst,Bst',zeros(2*nv,1);Bst,sparse(np,np),ones(np,1)/np; ...
            zeros(1,2*nv),ones(1,np)/np,zeros(1,1)]\[nlres;0];
   dxns=dxnsz(1:end-1); multiplier=dxnsz(end);
   elseif qmethod<=1,
%  dxns = -[Jnst,Bst';Bst,-nubeta*C]\nlres;
%----------------------------------- stabilized version
% %    dxnsz= -[Jnst,Bst',zeros(2*nv,1);Bst,-nubeta*C,ones(np,1)/np; ...
% %             zeros(1,2*nv),ones(1,np)/np,zeros(1,1)]\[nlres;0];
   dxnsz= -[Jnst,Bst',zeros(2*nv,1);Bst,-beta*muC,ones(np,1)/np; ...
            zeros(1,2*nv),ones(1,np)/np,zeros(1,1)]\[nlres;0];
   dxns=dxnsz(1:end-1); multiplier=dxnsz(end);
   end
   xns = flowsol + dxns;
% compute residual of new solution
   if     qmethod>1,  N = navier_q2(xy,mv,tpde,xns);
   elseif qmethod<=1, N = navier_q1(xy,ev,tpde,xns);
   end
%    Anst = tpvd*A + [N, sparse(nv,nv); sparse(nv,nv), N];
   Anst = A + [N, sparse(nv,nv); sparse(nv,nv), N];
   [Anst,Bst,fst,gst] = flowbc(Anst,B,f,g,xy,bound);
%    Bst = tpddi*Bst; gst = tpddi*gst;
   if     qmethod>1,  nlres = [Anst,Bst';Bst,sparse(np,np)]*xns-[fst;gst];
   elseif qmethod<=1, % % nlres = [Anst,Bst';Bst,-nubeta*C]*xns-[fst;gst];
       nlres = [Anst,Bst';Bst,-beta*muC]*xns-[fst;gst];
   end
   nlres_norm = norm(nlres);
   nnv=length(fst); soldiff=norm(xns(1:nnv)-flowsol(1:nnv));
   fprintf('nonlinear residual is %e',nlres_norm)
   fprintf('\n   velocity change is %e\n',soldiff)
   if testrun ~= 1
       % plot solution
%        flowplot(qmethod,xns,By,Bx,A,xy,xyp,x,y,bound,spc,66); drawnow;
       flowplot(qmethod,xns,By,Bx,L,xy,xyp,x,y,bound,spc,66); drawnow;
       %flowplot09(qmethod,xns,By,Bx,A,xy,xyp,x,y,bound,bndxy,bnde,obs,contourn,spc,166);drawnow
       pause(1)
   end
   flowsol = xns;
%% end of Newton iteration loop 
end
if nlres_norm <= nlres0_norm * tol_nl, 
   fprintf('\nfinished, nonlinear convergence test satisfied\n\n');
% explicitly reorthogonalize to remove hydrostatic pressure component   
   nlres = nlres - [zeros(2*nv,1);(sum(nlres(2*nv+1:2*nv+np))/np)*ones(np,1)];
else
   fprintf('\nfinished, stopped on iteration counts\n\n');
end
%
%%% estimate errors
if qmethod==1
   [jmpx,jmpy,els] = stressjmps_q1p0(viscosity,flowsol,xy,ev,ebound);
   [error_x,error_y,fex,fey,ae] = navierpost_q1p0_p(viscosity,flowsol,jmpx,jmpy,els,xy,ev);
   [error_x,error_y] = navierpost_q1p0_bc(viscosity,ae,fex,fey,...
                                        error_x,error_y,xy,ev,ebound);
   error_div = q1div(xy,ev,flowsol);
   errorest=sqrt(sum(error_x.^2 + error_y.^2 + error_div.^2));
   fprintf('estimated overall error is %10.6e \n',errorest)
   ee_error=sqrt((error_x.^2 + error_y.^2 + error_div.^2));
%% plot element errors
   %eplot(ee_error,ev,xy,x,y,67);
%% plot macroelement errors
   mplot(ee_error,ev,xy,x,y,67);  title('Estimated error')
   pause(5), figure(66)
elseif qmethod==0
   error_div = q1div(xy,ev,flowsol);	
elseif qmethod>1, %navierpost, 
end

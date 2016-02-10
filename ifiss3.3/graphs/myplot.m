existfignum = exist('fignum','var');
if ~existfignum
    fignum = 50;
else
    fignum = fignum+1;
end
whatplot = default('flowsol/x_it/flowsol_new 1/2/3 (default 3)',3);
if whatplot == 1
    flowplot(qmethod,flowsol,By,Bx,A,xy,xyp,x,y,bound,spc,fignum);
elseif whatplot == 2
    flowplot(qmethod,x_it,By,Bx,A,xy,xyp,x,y,bound,spc,fignum);
elseif whatplot == 3
    flowsol_new = flowsol-x_it;
    flowplot(qmethod,flowsol_new,By,Bx,A,xy,xyp,x,y,bound,spc,fignum);
else
    error('What to plot is not valid, entry should be 1 or 2.')
end
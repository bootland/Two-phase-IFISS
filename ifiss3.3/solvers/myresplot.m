function myresplot(resid)
%RESPLOT plot residuals computed by iterative solvers
%   resplot(resid);
%   input
%          resid        vector of residuals
%
%   IFISS function: AR, DJS, HCE, DJS; 17 January 2010.
% Copyright (c) 2005 D.J. Silvester, H.C. Elman, A. Ramage

load('testrun.mat')
% check current plot status
stest = get(gcf,'Children');
if size(stest,1)~=0,
   % ask for plotting information
   if testrun ~= 1
       figinfo=default('use new (enter figno) or existing (0) figure, default is 0',0);
   else
       if ( v == 1 && d == 1)
           figinfo = l + (r-1)*lL;
       else
           figinfo = 0;
       end
   end
   if figinfo==0,
       if testrun ~= 1
           figno=default('figure number (default is current active figure)',gcf);
       else
           figno = l + (r-1)*lL;
       end
      if figno==0, figno=gcf; end
      % re-use existing figure
      figure(figno); hold on
   elseif figinfo>0 && floor(figinfo)==figinfo,
      figure(figinfo); hold off;
   else
      fprintf('Illegal figure number, generating new figure.\n');
      hold off; figure;
   end
end

its=length(resid)-1;
if testrun ~= 1
    colval=default('colour (b,g,r,c,m,y,k,o): enter 1--8 (default 1)',1);
    lsval=default('linestyle (-,--,:,-.): enter 1--4 (default 1)',1);
else
    colval = v + (d-1)*lVISCD;
    if colval < 8.5
        lsval = 1;
    else
        lsval = 2;
        colval = colval + 2;
    end
end
Colour={'b','g','r','c','m',[1 0.84 0],'k',[0.85 0.33 0.1]};
remcol = rem(colval,length(Colour));
if remcol == 0
    remcol = length(Colour);
end
colour=Colour{remcol};
LS={'-','--',':','-.'};
remls = rem(lsval,length(LS));
if remls == 0
    remls = length(LS);
end
ls = LS{remls};
semilogy(0:its,resid,'color',colour,'LineStyle',ls,'LineWidth',2)
axis('square'), xlabel('Iterations'), ylabel(' log_{10}(residual)')
title('Residual reduction')
set(gca,'FontSize',14)
hold on
figure(gcf);

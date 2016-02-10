%navier_testproblem    sets up reference Examples 7.1 to 7.4
%   IFISS scriptfile: DJS; 29 September 2013.
% Copyright (c) 2009 D.J. Silvester, H.C. Elman, A. Ramage 
gohome
clear variables
load('testrun.mat')
if testrun == 1
    disp('Running test')
end
fprintf('\nspecification of reference Navier-Stokes problem.\n');
fprintf('\nchoose specific example (default is cavity)');
fprintf('\n     1  Channel domain');
fprintf('\n     2  Flow over a backward facing step');
fprintf('\n     3  Lid driven cavity');
fprintf('\n     4  Flow over a plate');
fprintf('\n     5  Flow over an obstacle\n');
if testrun ~= 1
    sn = default('',3);
else
    sn = 3;
end
if sn==1,
system('/bin/cp ./stokes_flow/test_problems/poiseuille_flow.m ./stokes_flow/specific_flow.m');
system('/bin/cp ./stokes_flow/test_problems/poiseuille_bc.m ./stokes_flow/stream_bc.m');
   channel_stokes; xsolve_navier;
elseif sn==2,
system('/bin/cp ./stokes_flow/test_problems/backwardstep_flow.m ./stokes_flow/specific_flow.m');
system('/bin/cp ./stokes_flow/test_problems/backwardstep_bc.m ./stokes_flow/stream_bc.m');
   longstep_stokes; solve_step_navier;
elseif sn==3,
    if testrun ~= 1
        lid_model=default('cavity type leaky/tight/regularised 1/2/3 (regularised)',3);
    else
        lid_model = 3;
    end
   if lid_model ==3,
system('/bin/cp ./stokes_flow/test_problems/regcavity_flow.m ./stokes_flow/specific_flow.m');
   elseif lid_model ==2,
system('/bin/cp ./stokes_flow/test_problems/tightcavity_flow.m ./stokes_flow/specific_flow.m');
   else
system('/bin/cp ./stokes_flow/test_problems/leakycavity_flow.m ./stokes_flow/specific_flow.m');
   end
system('/bin/cp ./stokes_flow/test_problems/zero_bc.m ./stokes_flow/stream_bc.m');
   if testrun ~= 1
       phases=default('one-phase/two-phase 1/2 (two-phase)',2);
   else
       phases = 2;
   end
   if phases==1
       square_stokes; solve_navier;
   elseif phases==2
       square_stokes_2; solve_navier_2;
   end
elseif sn==4,
system('/bin/cp ./stokes_flow/test_problems/plate_flow.m ./stokes_flow/specific_flow.m');
system('/bin/cp ./stokes_flow/test_problems/plate_bc.m ./stokes_flow/stream_bc.m');
   plate_stokes; solve_plate_navier;
elseif sn==5,
system('/bin/cp ./stokes_flow/test_problems/obstacle_flow.m ./stokes_flow/specific_flow.m');
system('/bin/cp ./stokes_flow/test_problems/obstacle_bc.m ./stokes_flow/stream_bc.m');
    obstacle_stokes; solve_obstacle_navier
else
   error('reference problem datafile not found!');
end

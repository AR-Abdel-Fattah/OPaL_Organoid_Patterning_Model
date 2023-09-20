%   OPaL_One
%   Please see README file First

size = 46;
nuc_wid = 6;
time = 10;
iterations = 500;
threshold = 2.36;
induction = 0.5;
ca = 0.050;
a  = 1;

%call OPaL_Single_Run Function
[Exp_scattering_organoids,Exp_patterning_organoids,Exp_negative_organoids,Table_of_outputs,Initial_sc_positions,Final_sc_positions,Initial_activation_signal,Final_activation_signal] = OPaL_Single_Run(size,nuc_wid,time,iterations,threshold,induction,ca,a);




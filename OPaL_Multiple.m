%   OPaL_Multiple
%   Please see README file First

size = 46;
nuc_wid = 6;
time = 10;
iterations = 500;
induction = 0.5;
a = 1;

threshold_ini = 0.0;
threshold_end = 5;
threshold_inc = 0.5;

ca_ini = 0.0;
ca_end = 0.5;
ca_inc = 0.25;


%call OPaL_Single_Run Function
[result_matrix_scattering,result_matrix_patterning,result_matrix_negative] = OPaL_Multiple_Run(size,nuc_wid,time,iterations,induction,a,threshold_ini,threshold_end,threshold_inc,ca_ini,ca_end,ca_inc);

 


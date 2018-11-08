function [n_amp_Sol, n_mero_Sol,n_mid_Sol,n_syn_Sol,n_mono_Sol, n_0_Sol] = kinetic_error_correction(k_on,k_off,nkt)
%{
k_on: MT capture rate
k_off: MT release rate 
nkt: number of KTs

n_amp: # ampheitelic chromosome
n_mero: # merotelic chromosome
n_mono: # monotellic chromosome
n_mid: # chromosomes in the middle of the graph
n_syn: # syntelic chromosome
n_0: # unattached chromosome
%}
syms n_amp(t) n_mero(t) n_mid(t) n_syn(t) n_mono(t) n_0(t)
ode1 = diff(n_amp) == k_on*n_mono + k_off*n_mero;
ode2 = diff(n_mero) == -3*k_off*n_mero+2*k_on*n_mid+2*k_on*n_syn;
ode3 = diff(n_mid) == -(2*k_off+2*k_on)*n_mid+k_on*n_mono+k_off*n_mero;
ode4 = diff(n_syn) == -(2*k_on+2*k_off)*n_syn+k_on*n_mono+k_off*n_mero;
ode5 = diff(n_mono) == -(3*k_on+k_off)*n_mono+(2*k_on*n_0+2*k_off*n_mid+2*k_off*n_syn)-k_on*n_amp;
ode6 = diff(n_0) == -2*k_on*n_0+k_off*n_mono;
odes=[ode1;ode2;ode3;ode4;ode5;ode6];

cond1 = n_0(0) == nkt;
cond2 = n_amp(0) == 0;
cond3 = n_mero(0) == 0;
cond4 = n_mid(0) == 0;
cond5 = n_syn(0) == 0;
cond6 = n_mono(0) == 0;
conds = [cond1; cond2; cond3; cond4; cond5; cond6];

[n_amp_Sol(t), n_mero_Sol(t),n_mid_Sol(t),n_syn_Sol(t),n_mono_Sol(t), n_0_Sol(t)] = dsolve(odes,conds);

end
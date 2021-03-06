function mpc = case30_E_Mishra
%CASE30    Power flow data for 30 bus, 6 generator case.
%   Please see CASEFORMAT for details on the case file format.
%
%   Based on data from ...
%     Alsac, O. & Stott, B., "Optimal Load Flow with Steady State Security",
%     IEEE Transactions on Power Apparatus and Systems, Vol. PAS 93, No. 3,
%     1974, pp. 745-751.
%   ... with branch parameters rounded to nearest 0.01, shunt values divided
%   by 100 and shunt on bus 10 moved to bus 5, load at bus 5 zeroed out.
%   Generator locations, costs and limits and bus areas were taken from ...
%     Ferrero, R.W., Shahidehpour, S.M., Ramesh, V.C., "Transaction analysis
%     in deregulated power systems using game theory", IEEE Transactions on
%     Power Systems, Vol. 12, No. 3, Aug 1997, pp. 1340-1347.
%   Generator Q limits were derived from Alsac & Stott, using their Pmax
%   capacities. V limits and line |S| limits taken from Alsac & Stott.

%   MATPOWER
%   $Id: case30.m 1559 2010-03-10 18:08:32Z ray $

%% MATPOWER Case Format : Version 2
mpc.version = '2';

%%-----  Power Flow Data  -----%%
%% system MVA base
mpc.baseMVA = 100;

%% bus data
% bus_i	type Pd     Qd      Gs      Bs	  area      Vm      Va    baseKV   zone Vmax	Vmin
mpc.bus = [
	1	3	0       0       0       0       1       1.060       0.0          135     1	1.05	0.95;
	2	2	21.7	12.7	0       0       1       1.043       -5.48        135     1	1.1     0.95;
	3	1	2.4     1.2     0       0       1       1.021       -7.96        135     1	1.05	0.95;
	4	1	7.6     1.6     0       0       1       1.012       -9.62        135     1	1.05	0.95;
	5	2	94.2    19      0       0.19	1       1.010       -14.37       135     1	1.05	0.95;
	6	1	0       0       0       0       1       1.010       -11.34       135     1	1.05	0.95;
	7	1	22.8	10.9	0       0       1       1.002       -13.12       135     1	1.05	0.95;
	8	2	30      30      0       0       1       1.010       -12.10       135     1	1.05	0.95;
	9	1	0       0       0       0       1       1.051       -14.38       135     1	1.05	0.95;
	10	1	5.8     2       0       0       3       1.045       -15.97       135     1	1.05	0.95;
	11	2	0       0       0       0       1       1.082       -14.39       135     1	1.05	0.95;
	12	1	11.2	7.5     0       0       2       1.057       -15.24       135     1	1.05	0.95;
	13	2	0       0       0       0       2       1.071       -15.24       135     1	1.1     0.95;
	14	1	6.2     1.6     0       0       2       1.042       -16.13       135     1	1.05	0.95;
	15	1	8.2     2.5     0       0       2       1.038       -16.22       135 	 1	1.05	0.95;
	16	1	3.5     1.8     0       0       2       1.045       -15.83       135     1	1.05	0.95;
	17	1	9       5.8     0       0       2       1.040       -16.14       135     1	1.05	0.95;
	18	1	3.2     0.9     0       0       2       1.028       -16.82       135     1	1.05	0.95;
	19	1	9.5     3.4     0       0       2       1.026       -17.00       135     1	1.05	0.95;
	20	1	2.2     0.7 	0       0       2       1.030       -16.80       135     1	1.05	0.95;
	21	1	17.5	11.2	0       0       3       1.033       -16.42       135 	 1	1.05	0.95;
	22	2	0       0       0       0       3       1.033       -16.41       135     1	1.1     0.95;
	23	1	3.2     1.6     0       0       2       1.027       -16.61   	 135     1	1.1     0.95;
	24	1	8.7     6.7 	0       0.043	3       1.021       -16.78       135     1	1.05	0.95;
	25	1	0       0       0       0       3       1.017       -16.35       135 	 1	1.05	0.95;
	26	1	3.5     2.3 	0       0       3       1.000       -16.77       135     1	1.05	0.95;
	27	1	0       0       0       0       3       1.023       -15.82       135     1	1.1     0.95;
	28	1	0       0       0       0       1       1.007       -11.97       135     1	1.05	0.95;
	29	1	2.4     0.9 	0       0       3       1.003       -17.06       135 	 1	1.05	0.95;
	30	1	10.6	1.9     0       0       3       0.992       -17.94       135     1	1.05	0.95;
];

%% generator data
%  bus	Pg      Qg      Qmax    Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
% DADOS ORIGINAIS
% mpc.gen = [
% 	1	260.2	-16.1	1000    -1000   1.060	100	1       200     50	0	0	0	0	0	0	0	0	0	0	0;
% 	2	40.0	50.0	50      -40     1.045	100	1       80      20	0	0	0	0	0	0	0	0	0	0	0;
% 	5	0.0 	37.0	40  	-40     1.010	100	1       50      15	0	0	0	0	0	0	0	0	0	0	0;
% 	8	0.0 	37.3	40  	-10     1.010	100	1       35      10	0	0	0	0	0	0	0	0	0	0	0;
% 	11	0.0 	16.2	24      -6      1.082	100	1       30      10	0	0	0	0	0	0	0	0	0	0	0;
% 	13	0.0     10.6	24  	-6      1.071	100	1       40      12	0	0	0	0	0	0	0	0	0	0	0;
% ];

% DADOS MISHRA 2015
mpc.gen = [                                                                                                         % TGER      vi   vr  v0      wr     d      Kp       Kr
	1	260.2	-16.1	250     -20     1.060	100	1       200     50	0	0	0	0	0	0	0	0	0	0	0   1            3   15  30      200    1      1       1;
	2	40.0	50.0	50      -40     1.045	100	1       80      20	0	0	0	0	0	0	0	0	0	0	0   1            3   15  30      80     1      1       1;
	5	0.0 	37.0	40  	-40     1.010	100	1       50      15	0	0	0	0	0	0	0	0	0	0	0   1            3   15  30      50     1      1       1;
	8	0.0 	37.3	40  	-10     1.010	100	1       35      10	0	0	0	0	0	0	0	0	0	0	0   1            3   15  30      35     1      1       1;
	11	0.0 	16.2	24      -6      1.082	100	1       30      10	0	0	0	0	0	0	0	0	0	0	0   1            3   15  30      30     1      1       1;
	13	0.0     10.6	24  	-6      1.071	100	1       40      12	0	0	0	0	0	0	0	0	0	0	0   1            3   15  30      40     1      1       1;
    22	0.0     10.6	24  	-20     1.071	100	1       40      0	0	0	0	0	0	0	0	0	0	0	0   2            3   10.28  25   40     0     2       10;
];


%% branch data
% fbus tbus	r       x       b	  rateA	rateB rateC	ratio	angle	status	angmin	angmax
mpc.branch = [
	1	2	0.0192	0.0575	0.0528	130	130	130     0       0	1	-360	360;
	1	3	0.0452	0.1652	0.0408	130	130	130     0       0	1	-360	360;
	2	4	0.0570	0.1737	0.0368	65	65	65      0       0	1	-360	360;
	3	4	0.0132	0.0379	0.0084  130	130	130     0       0	1	-360	360;
	2	5	0.0472	0.1983  0.0418	130	130	130     0       0	1	-360	360;
	2	6	0.0581	0.1763	0.0374	65	65	65      0       0	1	-360	360;
	4	6	0.0119	0.0414	0.0090  90	90	90      0       0	1	-360	360;
	5	7	0.0460	0.1160	0.0204	70	70	70      0       0	1	-360	360;
	6	7	0.0267	0.0820	0.0170	130	130	130     0       0	1	-360	360;
	6	8	0.0120	0.0420	0.0090  32	32	32      0       0	1	-360	360;
	6	9	0       0.2080	0       65	65	65      0.978   0	1	-360	360;
	6	10	0       0.5560	0       32	32	32      0.969   0	1	-360	360;
	9	11	0       0.2080	0       65	65	65      0       0	1	-360	360;
	9	10	0       0.1100	0       65	65	65      0       0	1	-360	360;
	4	12	0       0.2560	0       65	65	65      0.932   0	1	-360	360;
	12	13	0       0.1400	0       65	65	65      0       0	1	-360	360;
	12	14	0.1231	0.2559	0       32	32	32      0       0	1	-360	360;
	12	15	0.0662	0.1304	0       32	32	32      0       0	1	-360	360;
	12	16	0.0945	0.1987  0       32	32	32      0       0	1	-360	360;
	14	15	0.2210	0.1997  0       16	16	16      0       0	1	-360	360;
	16	17	0.0524	0.1923	0       16	16	16      0       0	1	-360	360;
	15	18	0.1073	0.2185	0       16	16	16      0       0	1	-360	360;
	18	19	0.0639	0.1292	0       16	16	16      0       0	1	-360	360;
	19	20	0.0340	0.0680	0       32	32	32      0       0	1	-360	360;
	10	20	0.0936	0.2090	0       32	32	32      0       0	1	-360	360;
	10	17	0.0324	0.0845	0       32	32	32      0       0	1	-360	360;
	10	21	0.0348	0.0749	0       32	32	32      0       0	1	-360	360;
	10	22	0.0727	0.1499	0       32	32	32      0       0	1	-360	360;
	21	22	0.0116	0.0236	0       32	32	32      0       0	1	-360	360;
	15	23	0.1000  0.2020  0       16	16	16      0       0	1	-360	360;
	22	24	0.1150	0.1790	0       16	16	16      0       0	1	-360	360;
	23	24	0.1320	0.2700	0       16	16	16      0       0	1	-360	360;
	24	25	0.1885	0.3292	0       16	16	16      0       0	1	-360	360;
	25	26	0.2544	0.3800	0       16	16	16      0       0	1	-360	360;
	25	27	0.1093	0.2087	0       16	16	16      0       0	1	-360	360;
	28	27	0       0.3960  0       65	65	65      0.968   0	1	-360	360;
	27	29	0.2198	0.4153	0       16	16	16      0       0	1	-360	360;
	27	30	0.3202	0.6027  0       16	16	16      0       0	1	-360	360;
	29	30	0.2399	0.4533	0       16	16	16      0       0	1	-360	360;
	8	28	0.0636	0.2000  0.0428	32	32	32      0       0	1	-360	360;
	6	28	0.0169	0.0599	0.0130	32	32	32      0       0	1	-360	360;
];

%%-----  OPF Data  -----%%
%% area data
%	area	refbus
mpc.areas = [
	1	8;
	2	23;
	3	26;
];

%% generator cost data
%	1	startup	shutdown	n	x1	y1	...	xn	yn
%	2	startup	shutdown	n	c(n-1)	...	c0

% DADOS ORIGINAIS
% mpc.gencost = [
% 	2	0	0	3	0.02	2       0;
% 	2	0	0	3	0.0175	1.75	0;
% 	2	0	0	3	0.0625	1       0;
% 	2	0	0	3	0.00834	3.25	0;
% 	2	0	0	3	0.025	3       0;
% 	2	0	0	3	0.025	3       0;
% ];

% % DADOS BASEADOS NA DISSERTA??O DE MESTRADO DA ELIS PG 136 (160 DO PDF)
% mpc.gencost = [     % a($/MW^2)  b($/MW) c($)   e       f
% 	2	0	0	3	0.00375     2       0       120     0.073;  
% 	2	0	0	3	0.0175      1.75	0       50      0.032;
% 	2	0	0	3	0.0625      1       0       30      0.051;
% 	2	0	0	3	0.00834     3.25	0       25      0.026;
% 	2	0	0	3	0.025       3       0       25      0.026;
% 	2	0	0	3	0.025       3       0       30      0.048;
% ];

% DADOS BASEADOS NO ARTIGO RICARDO 2019
% mpc.gencost = [     % a($/MW^2)  b($/MW) c($)   e       f
% 	2	0	0	3	0.00375     2       0       18     0.037;  
% 	2	0	0	3	0.0175      1.75	0       16      0.038;
% 	2	0	0	3	0.0625      1       0       14      0.04;
% 	2	0	0	3	0.00834     3.25	0       12      0.045;
% 	2	0	0	3	0.025       3       0       13      0.042;
% 	2	0	0	3	0.025       3       0       13.5      0.041;
%     2	0	0	3	0.025       3       0       13      0.042
% ];

% DADOS BASEADOS NO ARTIGO DE MISHRA 2015 
mpc.gencost = [     % a($/MW^2)  b($/MW) c($)   e       f
	2	0	0	3	0.00375     2       0       0     0;  
	2	0	0	3	0.0175      1.75	0       0      0;
	2	0	0	3	0.0625      1       0       0      0;
	2	0	0	3	0.00834     3.25	0       0      0;
	2	0	0	3	0.025       3       0       0      0;
	2	0	0	3	0.025       3       0       0      0;
    2	0	0	3	0.025       3       0       0      0
];


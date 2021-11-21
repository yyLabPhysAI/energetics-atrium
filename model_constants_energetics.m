data.Is 		= -15000;           % [???]			STIMULUS
data.R     		= 8314;          	% [mJ/(mol*K)]	Ideal gas const.      
data.T     		= 308;           	% [K]			Temp.                 
data.Cm    		= 0.05;          	% [pF] 			Membrane capacitance  
data.F     		= 96487;         	% [C/mol]		Faraday const         
data.RTONF 		= data.T*0.08554;	% [???]			RT/F 					
data.Ca_i   	= 0.000071;     	% [???]			Ca^+2 	in 	        
data.Ca_o  		= 2.5;           	% [???]			Ca^+2 	out         
data.Na_i   	= 8.4;          	% [mM]			Na^+ 	in 	        
data.Na_o   	= 140;				% [mM]			Na^+	out	        
data.K_i    	= 140;				% [mM]			K^+		in	        
data.K_o    	= 5.0;				% [mM]			K^+		out	        
data.V_i    	= 0.0126;			% [nL]			Intracellular volume
data.V_Ca   	= 0.005884;			% [nL]			Ca2+ intracellular volume
data.V_c    	= 0.0025;			% [nL]			Cytoplasmic volume
data.V_up   	= 0.0003969;		% [nL]			SR uptake compartment volume
data.V_rel  	= 0.000044;			% [nL]			SR release compartment volume
data.V_nSR  	= 0.0408;			% [nL]			???			        

%%Iup
data.K_cyCa		= 0.0003*10^3;		% []			Equilibrium binding ca2+
data.K_xcs		= 0.4;				% []			Translocation constant
data.K_srCa		= 0.5*10;			% []			Equilibrium binding ca2+ concentration on uptake of compartment of the SR side

% force parameters
data.F_f    	= 0.04*1e3;			% [1/ms] 		???, corrected from msec to sec ??? so is it ms or sec???
data.SL_0  		= 0.8;				% [um] 			Minimal length of the sarcomere
data.N_c   		= 2e13;				% [1/mm^2] 		The reciprocal of the crosssection of the tissue
data.F_k0  		= 350;				% [1/mM] 		The crossbridge independent coefficent of calcium affinity
data.F_k1  		= 3000;				% [1/mM] 		Cooperativity coefficient 
data.F_k05 		= 2.5e9;			% [1/mm^3]		The half-maximal crossbridge calcium affinity
data.F_kl  		= 60*1e3;			% [1/mM/ms] 	Rate constant of calcium binding to troponin low-affinity sites, corrected from msec to sec ??? so is it ms or sec???
data.F_go  		= 0.03*1e3;			% [1/ms]		The crossbridge weakening rate at isometric regime, corrected from msec to sec ??? so is it ms or sec???
data.F_gl  		= 4.4e6; 			% m[1/m] 		Mechanical feedback coefficient
data.F_XB  		= 2e-9;				% [mN] 			Unitary force per crossbridge at isometric regime
data.FN 		= 3.5;				% [] 			Hill's coefficient

% energy consumption parameters
data.C_ATPi		= 2.6;				% [mM]			ATP consumption in the cytoplasm
data.MaxATP 	= 0.02533; 			% [mM]			Maximal ATP consumption by the sarcomeres
data.K_ACI  	= 0.016;			% [1/min]		Non-ca2+ AC activity
data.K_AC   	= 0.0735; 			% [1/min]		Non-ca2+ AC activation
data.K_Ca   	= 0.000178;			% [mM]			Maximal Ca2+ AC activation
data.k_bCM  	= 0.5421;			% [1/ms]		Ca2+ dissociation constant for calmodulin
data.k_fCM  	= 227.7;			% [1/(mM*ms)]	Ca2+ association constant for calmodulin
data.K_ACCa 	= 0.000024;			% [M]			Half-maximal ca2+ AC activation
data.k_XBATP	= 0.03;				% []			Coefficient of force generation at the energy state of the cell ??? why the same comment as k_XBADP ???
data.k_XBADP	= 0.26;				% []			Coefficient of force generation at the energy state of the cell 

% Mitochondrial calcium handling parameters
data.V_uni_max	= 0.0275;			% [mM/ms] 		Vmax uniporter Ca2+ transport
data.psi0		= 91; 				% [mV]			offset membrane potential
data.K_act		= 3.8e-4; 			% [mM]			activation constant
data.K_trans	= 0.019;			% [mM]			Kd for translocated Ca2+
data.L			= 110; 				% [] 			Keq for conformational transitions in uniporter
data.n_a		= 2.8;			 	% [] 			uniporter activation cooperativity
data.V_NaCa_max	= 10^-4;			% [mM/ms]		Vmax of Na/Ca antiporter
data.b			= 0.5;			 	% [] 			deltaPsi_m dependence of Na/Ca antiporter 
data.K_Na		= 8; 				% [mM]			antiporter Na+ constant
data.K_Ca		= 0.008;			% [mM]			antiporter Ca2+ constant
data.n			= 3; 				% []			Na/Ca2+ antiporter cooperativity
data.delta		= 3e-4;			 	% [] 			Fraction of free [Ca2+]m
 
% Oxidative phosphorylation parameters
data.r_a		= 6.394e-13;		% 1/[ms] 		sum of products of rate constants
data.r_c1		= 2.656e-22*10^5; 	% 1/[ms] 		sum of products of rate constants
data.r_c2		= 8.632e-30;		% 1/[ms] 		sum of products of rate constants
data.r_1		= 2.077e-18*5*10^4/2.8415; % [???] 	sum of products of rate constants
data.r_2		= 1.728e-9;		 	% [???]			sum of products of rate constants
data.r_3		= 1.059e-26;		% [???]			sum of products of rate constants
data.rho_res	= 3e-3*10^5; 		% [mM] 			concentration of electron carriers (respiratory complexes I-III-IV)
data.K_res		= 1.35e18;			% []   			equilibrium constant of respiration
data.rho_resF	= 3.75e-4; 			% [mM] 			concentration of electron carriers (respiratory complexes II-III-IV)
data.psi_B		= 50*30; 			% [mV] 			phase boundary potential
data.G			= 0.85; 			% [] 			correction factor for voltage
data.K_resF		= 5.765e13; 		% [] 			equilibrium constant FADH2 oxidation
data.r_b		= 1.76e-16; 		% 1/[ms]		sum of products of rate constants
data.FADH2		= 1.24; 			% [mM] 			concentration of FADH2 (reduced)
data.FADH		= 0.02; 			% [mM] 			concentration of FAD (oxidized)
data.p_a		= 1.656e-8;			% 1/[ms]		sum of products of rate constants
data.p_b		= 3.3373e-10;	 	% 1/[ms]		sum of products of rate constants
data.p_c1		= 9.651e-17*10^10;	% 1/[ms]		sum of products of rate constants
data.p_c2		= 7.739e-17*10^10;	% 1/[ms]		sum of products of rate constants
data.p_1		= 1.346e-8;	 		% [] 			sum of products of rate constants
data.p_2		= 7.739e-7;	 		% [] 			sum of products of rate constants
data.p_3		= 6.65e-15*10^8; 	% [] 			sum of products of rate constants
data.KCaATP		= 165/10^6;			% [???]			???
data.rho_F1		= 1.5;				% [mM] 			concentration of F1F0-ATPASE
data.K_F1		= 1.71e6;			% [] 			equilibrium constant of ATP hydrolysis
data.Pi			= 2; 				% [mM] 			inorganic phosphate concentration
data.C_A		= 1.5; 				% [mM] 			total sum of mitochondrial adenine nucleotides
data.V_ant_max	= 0.025/1400;		% [mM/ms] 		maximal rate of ANT
data.h_ANT		= 0.5/10; 			% [] 			fraction of delta_psi_m
data.g_H		= 10^-8		; 		% [mM]/([ms][mV])	ionic conductance of the inner membrane
data.delta_pH	= -0.6		; 		% [pH units] 	pH gradient across the inner membrane
data.C_PN		= 10;		 		% [mM] 			total sum of mitochondrial pyridine nucleotides
data.C_mito		= 1.812e-3; 		% [mM]/[mV] 	inner membrane capacitance
data.Cc			= 25;				% [???]			???

% Tricarbolxylic acid cycle parameters
data.C_AcCoA	= 1; 				% [mM]			acetyl coA concentration
data.k_cat_cs	= 0.05;				% [1/ms] 		catalytic constant of CS
data.E_T_cs		= 0.4;				% [mM]			concentration of CS
data.K_M_AcCoA	= 1.26e-2*100; 		% [mM]			michaelis constant of AcCoA
data.K_M_OAA	= 6.4e-4;			% [mM]			michaelis constant of OAA
data.C_K_int	= 1; 				% [mM]			Sum of TCA cycle intermediates' concentration
data.k_f_ACO	= 1.25e-2; 			% 1/[ms]		forward rate constant of ACO
data.k_E_ACO	= 2.22/1000; 		% [] 			equilibrium constant of ACO
data.K_ADP_a	= 0.62;				% [mM]			activation constant by ADP
data.K_Ca_a		= 0.0005;			% [mM]			activation constant by ADP
data.K_i_NADH	= 0.19;				% [mM]			inhibition constant by NADH
data.k_cat_IDH	= 0.03;				% 1/[ms]		rate constant of IDH
data.E_T_IDH	= 0.109;			% [mM]			concentration of IDH
data.C_H		= 2.5e-5; 			% [mM]			matrix proton concentration
data.k_h_1		= 8.1e-5; 			% [mM]			ionization constant of IDH
data.k_h_2		= 5.98e-5;			% [mM]			ionization constant of IDH
data.K_M_ISOC	= 1.52;				% [mM]			michaelis constant for isocitrate
data.N_i		= 2; 				% [] 			cooperativity of isocitrate
data.K_M_NAD	= 0.923;			% [mM]			michaelis constant for NAD+
data.K_M_Mg		= 0.0308;			% [mM]			activation constant for Mg2+
data.K_M_Ca		= 1.27e3; 			% [mM]			activation constant for Ca2+
data.E_T_KGDH	= 0.5; 				% [mM]			concentration of KGDH
data.k_cat_KGDH	= 0.05; 			% 1/[ms]		rate constant of KGDH
data.K_M_aKG	= 1.94/1000;		% [mM]			michaelis constant for aKG
data.K_M_NAD_new= 38.7/1000;		% [mM]			michaelis constant for NAD
data.n_aKG		= 1.2/2.5; 			% [] 			Hill coefficient of KGDH for aKG
data.C_Mg		= 0.4;  			% [mM]			Mg2+ concentration in mitochondria
data.k_f_SL		= 5e-4*10^6; 		% [mM]/[ms]		forward rate constant of SL
data.k_E_SL		= 3.115;  			% [] 			equilibrium constant of the SL reaction
data.C_CoA		= 0.02/100; 		% [mM]			coenzyme A cocentration
data.k_cat_SDH	= 3e-3;  			% 1/[ms]		rate constant of SDH
data.E_T_SDH	= 0.5;  			% [mM]			SDH enzyme concentration
data.K_M_SUC	= 0.03/1000; 		% [mM]			michaelis constant for succinate
data.K_i_FUM	= 1.3;  			% [mM]			inhibition constant by oxalacetate
data.K_i_sdh_OAA= 0.15;  			% [mM]			inhibition constant by oxalacetate
data.k_f_FH		= 3.32e-3; 			% 1/[ms]		forward rate constant of FH
data.K_E_FH		= 1/1000;  			% [] 			equilibrium constant of FH
data.k_h1		= 1.13e-5; 			% [mM]			ionization constant of MDH
data.k_h2		= 26.7;  			% [mM]			ionization constant of MDH
data.k_h3		= 6.68e-9; 			% [mM]			ionization constant of MDH
data.k_h4		= 5.62e-6; 			% [mM]			ionization constant of MDH
data.k_offset	= 3.99e-2; 			% [] pH			independent term in the pH activation factor of MDH
data.k_cat_MDH	= 0.111;  			% 1/[ms] 		rate constant of MDH
data.E_T_MDH	= 0.154;  			% [mM]			Total MDH enzyme concentration
data.K_M_MAL	= 1.493;  			% [mM]			michaelis constant for malate
data.K_i_OAA	= 3.1e-3; 			% [mM]			inhibition constant for oxalacetate
data.K_M_NAD_mdh= 0.2244;  			% [mM]			michaelis constant for NAD+
data.C_GLU		= 10; 				% [mM]			glutamate concentration
data.k_f_AAT	= 6.44e-4*10000; 	% 1/[ms]		forward rate constant of AAT
data.K_E_AAT	= 6.6*10; 			% []			eqilibrium constant of AAT
data.k_ASP		= 1.5e-6*10^6; 		% 1/[ms]		rate constant of aspartate consumption
 
%cytoplasmic energy handling parameters
data.C_T		= 25;				% [mM]			total concentration of creatine metabolite (both compartments)
data.k_CK_cyto	= 1.4e-4;			% 1/[ms]		forward rate constant of cytoplasmic CK
data.k_CK_mito	= 1.33e-6;			% 1/[ms]		forward rate constant of mitochondrial CK
data.k_tr_Cr	= 2e-3;				% 1/[ms]		transfer rate constant of CrP
data.K_EQ		= 0.0095;			% [] 			equilibrium constant of CK
data.V_ATPase_cyto	= 10^-5;		% [???]			constitutive cytosolic ATP consumption rate
 
data.V_myo		= 25.84; 			% [pL]			cytosolic volume
data.V_mito		= 15.89; 			% [pL]			mitochondrial volume
data.A_cap		= 1.534e-4; 		% [cm^2] 		capacitative cell surface area

%%mitochondrial calcium
data.P_Ca 		= 2.159;			% 1/[ms]		permeability of calcium uniporter
data.Z_Ca 		= 2; 				% []			calcium valence
data.alpha_m 	= 0.2;				% []			mitochondrial calcium activity coefficient
data.alpha_e 	= 0.341;			% []			extramitochondrial calcium activity coefficient
data.V_NC 		= 1.863*10^-2;		% [mM/ms]			NCX maximal velocity

data.Na_e 		= 5;				% [mM]			extramitochdrial na+ concentration
data.Na_m 		= 3.96;				% [mM]			mitochondrial na+ concentration
data.beta_Ca 	= 0.1;				% []			Ca2+ fraction that binds to Ca2+ buffers in the mitochondria
data.K_D_Ca 	= 1.27*10^-3;		% [mM]			KGDHC Ca2+ binding constant
data.K_D_Mg 	= 0.0308;			% [mM]			KGDHC Mg2+ binding constant


%%force generation model
data.SL_o		= 0.8;				% [um]			minimal length of the sarcomere
data.N_c		= 2*10^13;			% [1/mm^2]		the reciprocal of the crosssection of the tissue
data.F_k0		= 350;				% [1/mM]		the crossbridge independent coefficent of calcium affinity
data.F_k1		= 3000;				% [1/mM]		cooperativity coefficient 
data.FN			= 3.5;				% []			Hill's coefficient
data.F_k_half	= 2.5*10^9;			% [1/mm^3]		the half-maximal crossbridge calcium affinity
data.F_kl		= 60;				% [1/mM/ms]		corrected from msec to sec %rate constant of calcium binding to troponin low-affinity sites
data.F_f		= 0.04;				% [1/ms]		corrected from msec to sec
data.F_g0		= 0.03;				% [1/ms]		corrected from msec to sec %the crossbridge weakening rate at isometric regime
data.F_gl		= 4.4*10^6;			% m[1/m]		mechanical feedback coefficient

data.K_M_ATP	= 0.03;				% []			coefficient of force generation at the energy state of the cell
data.K_M_ADP	= 0.26;				% []			coefficient of force generation at the energy state of the cell
data.Max_ATP	= 12e-4;			% [mM]			maximal ATP consumption by the sarcomeres
data.CATPi		= 2.6;				% [mM]			ATP consumption in the cytoplasm

%% Initial conditions for the state vector
y0(1) 	= -80;     					% V 			[mV] 		Transmembrane potential 
y0(2) 	= 0.00016; 					% pa 			[]			activation gating parameter for IK,r
y0(3) 	= 0.76898; 					% pi 			[]			inactivation gating parameter for IK,r
y0(4) 	= 0.02032; 					% n 			[]			activation gating parameter for IK,s
y0(5) 	= 0.00006; 					% r1 			[]			activation gating variable of It
y0(6) 	= 0.5753;  					% s1 			[]			Fast inactivation gating variable for It
y0(7) 	= 0.39871; 					% s2 			[]			Slow inactivation gating variable for It
y0(8) 	= 0.57363; 					% s3 			[]			Third inactivation gating variable for It
y0(9) 	= 0.01309; 					% m 			[]			activation gating variable for INa
y0(10)	= 0.706;   					% h1 			[]			Fast inactivation gating variable for INa
y0(11)	= 0.61493; 					% h2 			[]			Slow inactivation gating variable for INa
y0(12)	= 0.00003; 					% dL 			[]			activation gating variable for ICaL
y0(13)	= 0.99981; 					% fL 			[]			inactivation gating variable for ICaL
y0(14)	= 0.00046; 					% dT 			[]			activation gating variable for ICaT
y0(15)	= 0.30752; 					% fT 			[]			inactivation gating variable for ICaT

% state variables for LINDBLAD Ca++ handling
y0(16) 	= 8.4;       				% Na_i			[mM]		Intracellular Na+ concentration	
y0(17) 	= 0.730866;  				% Ca_up			[mM]		Ca2+ concentration in uptake compartment	
y0(18) 	= 0.726776;  				% Ca_rel		[mM]		Ca2+ concentration in release compartment	
y0(19) 	= 0.00007305;  				% Cai			[mM]		intracelluar Ca2+ concentration	
y0(20) 	= 0.029108;  				% fac			[]			fractional occupancy of calmodulin by Ca2+	
y0(21) 	= 0.014071;  				% faTc			[]			fractional occupancy of troponin-Ca2+ complex by ca2+	
y0(22) 	= 0.214036;  				% faTmgc		[]			fractional occupancy of troponin-Mg2+ complex by ca2+	
y0(23) 	= 0.693565;  				% faTmgm		[]			fractional occupancy of troponin-Mg2+ complex by Mg2+	
y0(24) 	= 0.465921;  				% faCalse		[]			fractional occupancy of caldsequestrin by Ca2+	
y0(25) 	= 5.0;       				% Kc			[mM]		extracellular K+ concentration	
y0(26) 	= 140.0;   					% Ki			[mM]		intracellular K+ concentration	// 100.0 - Lindblad ???? relevent???
y0(27) 	= 0.288039;					% F1			[]			Relative amount of inactive precursor in SR release compartment
y0(28) 	= 0.002262;					% F2			[]			Relative amount of activator in SR release compartment
y0(29) 	= 0.612697;					% F3			[]			Relative amount of inactive product in SR release compartment
y0(30) 	= 0.7755;  					% fca			[???]		???		<------ assumed constant in this use of the model, implemented for future use	<------ assumed constant in this use of the model, implemented for future use

% state variables for the force generation model
y0(31) 	= 1.75;    					% SL  			[um]		sarcomere length
y0(32) 	= 0.06;    					% A  			[]			density of regulatory units with bound ca2+ and adjacent weak crossbridges
y0(33) 	= 0.02;    					% TT  			[]			density of regulatory units with bound ca2+ and adjacent strong crossbridges
y0(34) 	= 0.06;    					% U  			[]			density of regulatory units without bound ca2+ but with adjacent strong crossbridges
y0(35) 	= 0;       					% Ve  			[um/sec]	velocity <------ assumed constant in this use of the model, implemented for future use

% state variables for the energy consumption model
y0(36) 	= 7.977;       				% ATPi 			[mM]		EC-coupling linked ATP concentration <------ assumed constant in this use of the model, implemented for future use

% state variable for the Calmodulin binding, needed for modeling of
% Adenylate cyclase
% y0(37) = 0.042;
y0(37)	= 0.000076; 				% Ca2+m 		[microM] 	mitochondrial free Ca2+ concentration 
y0(38)	= 7.972; 					% [ATP]ic 		[mM]		cytosolic ATP concentration not linked to EC coupling
y0(39)	= 18.297; 					% [CrP]i 		[mM]		mitochondrial-linked creatine phosphate concentration
y0(40)	= 18.291; 					% [CrP]ic 		[mM]		cytosolic creatine phosphate concentration
y0(41)	= 0.276; 					% [ADP]m 		[mM]		mitochondrial ADP concentration
y0(42)	= 5.403; 					% [NADH] 		[mM]		mitochondrial NADH concentration
y0(43)	= -140.7; 					% [deltaPsi_m] 	[mV]		Inner mitochondrial membrane potential
y0(44)	= 0.41; 					% [ISOC]		[mM]		Isocitrate concentration (mitochondrial)
y0(45)	= 2.596e-4;					% [aKG]			[mM]		Α-ketoglutarate concentration (mitochondrial)
y0(46)	= 0.362; 					% [SCoA]  		[mM]		Succinyl-CoA concentration (mitochondrial)
y0(47)	= 1.06e-4;					% [Suc]  		[mM]		Succinate concentration (mitochonrial)
y0(48)	= 0.0282; 					% [FUM] 		[mM]		fumarate concentration (mitochondrial)
y0(49)	= 0.01316; 					% [MAL] 		[mM]		Malate concentration (mitochondrial)
y0(50)	= 1.623e-6;					% [OAA] 		[mM]		Oxaloacetate concentration (mitochondrial)
y0(51)	= 5.403; 					% [FLV] 		[mM]		Flavoprotein concentration
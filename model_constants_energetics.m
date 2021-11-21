
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
data.V_i    	= 0.0126;			% [???]			Intracellular volume	
data.V_Ca   	= 0.005884;			% [mV]			???     		
data.V_c    	= 0.0025;			% [???]			Cytoplasmic volume	
data.V_up   	= 0.0003969;		% [mV]			???			        
data.V_rel  	= 0.000044;			% [mV]			???			        
data.V_nSR  	= 0.0408;			% [mV]			???			        

%%Iup
data.K_cyCa		= 0.0003*10^3;		% [???]			???
data.K_xcs		= 0.4;				% [???]			???
data.K_srCa		= 0.5*10;			% [???]			???

% force parameters
data.F_f    	= 0.04*1e3;			% [???]			???
data.SL_0  		= 0.8;				% [???]			???
data.N_c   		= 2e13;				% [???]			???
data.F_k0  		= 350;				% [???]			???
data.F_k1  		= 3000;				% [???]			???
data.F_k05 		= 2.5e9;			% [???]			???
data.F_kl  		= 60*1e3;			% [???]			???
data.F_go  		= 0.03*1e3;			% [???]			???
data.F_gl  		= 4.4e6; 			% [???]			???
data.F_XB  		= 2e-9;				% [???]			???
data.FN 		= 3.5;				% [???]			???

% energy consumption parameters
data.C_ATPi		= 2.6;				% [???]			???
data.MaxATP 	= 0.02533; 			% [???]			???
data.K_ACI  	= 0.016;			% [???]			???
data.K_AC   	= 0.0735; 			% [???]			???
data.K_Ca   	= 0.000178;			% [???]			???
data.k_bCM  	= 0.5421;			% [???]			???
data.k_fCM  	= 227.7;			% [???]			???
data.K_ACCa 	= 0.000024;			% [???]			???
data.k_XBATP	= 0.03;				% [???]			???
data.k_XBADP	= 0.26;				% [???]			???

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
data.P_Ca 		= 2.159;			%1/[ms]
data.Z_Ca 		= 2; 				% [???]			???
data.alpha_m 	= 0.2;				% [???]			???
data.alpha_e 	= 0.341;			% [???]			???
data.V_NC 		= 1.863*10^-2;		% [???]			???

data.Na_e 		= 5;				% [???]			???
data.Na_m 		= 3.96;				% [???]			???
data.beta_Ca 	= 0.1;				% [???]			???
data.K_D_Ca 	= 1.27*10^-3;		% [???]			???
data.K_D_Mg 	= 0.0308;			% [???]			???


%%force generation model
data.SL_o		= 0.8;				% [???]			???
data.N_c		= 2*10^13;			% [???]			???
data.F_k0		= 350;				% [???]			???
data.F_k1		= 3000;				% [???]			???
data.FN			= 3.5;				% [???]			???
data.F_k_half	= 2.5*10^9;			% [???]			???
data.F_kl		= 60;				% [???]			???
data.F_f		= 0.04;				% [???]			???
data.F_g0		= 0.03;				% [???]			???
data.F_gl		= 4.4*10^6;			% [???]			???
data.K_M_ATP	= 0.03;				% [???]			???
data.K_M_ADP	= 0.26;				% [???]			???
data.Max_ATP	= 12e-4;			% [???]			???
data.CATPi		= 2.6;				% [???]			???

%% Initial conditions for the state vector

y0(1) 	= -80;     					% V 			[mV] 		Transmembrane potential 
y0(2) 	= 0.00016; 					% pa 			[???]		???
y0(3) 	= 0.76898; 					% pi 			[???]		???
y0(4) 	= 0.02032; 					% n 			[???]		???
y0(5) 	= 0.00006; 					% r1 			[???]		???
y0(6) 	= 0.5753;  					% s1 			[???]		???
y0(7) 	= 0.39871; 					% s2 			[???]		???
y0(8) 	= 0.57363; 					% s3 			[???]		???
y0(9) 	= 0.01309; 					% m 			[???]		???
y0(10)	= 0.706;   					% h1 			[???]		???
y0(11)	= 0.61493; 					% h2 			[???]		???
y0(12)	= 0.00003; 					% dL 			[???]		???
y0(13)	= 0.99981; 					% fL 			[???]		???
y0(14)	= 0.00046; 					% dT 			[???]		???
y0(15)	= 0.30752; 					% fT 			[???]		???

% state variables for LINDBLAD Ca++ handling
y0(16) 	= 8.4;       				% Na_i			[???]		???	
y0(17) 	= 0.730866;  				% Ca_up			[???]		???	
y0(18) 	= 0.726776;  				% Ca_rel		[???]		???	
y0(19) 	= 0.00007305;  				% Cai			[???]		???	
y0(20) 	= 0.029108;  				% fac			[???]		???	
y0(21) 	= 0.014071;  				% faTc			[???]		???	
y0(22) 	= 0.214036;  				% faTmgc		[???]		???	
y0(23) 	= 0.693565;  				% faTmgm		[???]		???	
y0(24) 	= 0.465921;  				% faCalse		[???]		???	
y0(25) 	= 5.0;       				% Kc			[???]		???	
y0(26) 	= 140.0;   					% Ki			[???]		???	// 100.0 - Lindblad ???? relevent???
y0(27) 	= 0.288039;					% F1			[???]		???
y0(28) 	= 0.002262;					% F2			[???]		???
y0(29) 	= 0.612697;					% F3			[???]		???
y0(30) 	= 0.7755;  					% fca			[???]		???	<------ assumed constant in this use of the model, implemented for future use

% state variables for the force generation model
y0(31) 	= 1.75;    					% SL  			[???]		???
y0(32) 	= 0.06;    					% A  			[???]		???
y0(33) 	= 0.02;    					% TT  			[???]		???
y0(34) 	= 0.06;    					% U  			[???]		???
y0(35) 	= 0;       					% Ve  			[???]		??? <------ assumed constant in this use of the model, implemented for future use

% state variables for the energy consumption model
y0(36) 	= 7.977;       				% ATPi 			[???]		??? <------ assumed constant in this use of the model, implemented for future use

% state variable for the Calmodulin binding, needed for modeling of
% Adenylate cyclase

% y0(37) = 0.042;
y0(37)	= 0.000076; 				% Ca2+m 		[microM] 	???
y0(38)	= 7.972; 					% [ATP]ic 		[???]		???
y0(39)	= 18.297; 					% [CrP]i 		[???]		???
y0(40)	= 18.291; 					% [CrP]ic 		[???]		???
y0(41)	= 0.276; 					% [ADP]m 		[???]		???
y0(42)	= 5.403; 					% [NADH] 		[???]		???
y0(43)	= -140.7; 					% [deltaPsi_m] 	[???]		???
y0(44)	= 0.41; 					% [ISOC]		[???]		???
y0(45)	= 2.596e-4;					% [aKG]			[???]		???
y0(46)	= 0.362; 					% [SCoA]  		[???]		???
y0(47)	= 1.06e-4;					% [Suc]  		[???]		???
y0(48)	= 0.0282; 					% [FUM] 		[???]		???
y0(49)	= 0.01316; 					% [MAL] 		[???]		???
y0(50)	= 1.623e-6;					% [OAA] 		[???]		???
y0(51)	= 5.403; 					% [FLV] 		[???]		???

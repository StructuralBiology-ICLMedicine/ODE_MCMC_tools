#ifdef SIZE_DEFINITIONS
#define N_METABS 29
#define N_ODE_METABS 0
#define N_INDEP_METABS 22
#define N_COMPARTMENTS 1
#define N_GLOBAL_PARAMS 38
#define N_KIN_PARAMS 0
#define N_REACTIONS 34

#define N_ARRAY_SIZE_P  39	// number of parameters
#define N_ARRAY_SIZE_X  22	// number of initials
#define N_ARRAY_SIZE_Y  16	//14	//13	//11	//9	//7	// number of assigned elements
#define N_ARRAY_SIZE_XC 22	// number of x concentration
#define N_ARRAY_SIZE_PC 0	// number of p concentration
#define N_ARRAY_SIZE_YC 16	//14	//13	//11	//9	//7	// number of y concentration
#define N_ARRAY_SIZE_DX 22	// number of ODEs 
#define N_ARRAY_SIZE_CT 16	//14	//13	//11	//9	//7	// number of conserved totals

#endif // SIZE_DEFINITIONS

#ifdef TIME
#define T  <set here a user name for the time variable> 
#endif // TIME

#ifdef NAME_ARRAYS
const char* p_names[] = {"TestTube", "tl_rate_gfp", "tl_rate_rfp", "tl_escape_rate_gfp", "tl_escape_rate_rfp", "tx_rate", "tx_escape_rate_gfp", "tx_escape_rate_rfp", "rna_deg_rate", "gfp_mat_rate", "rfp_mat_rate1", "rfp_mat_rate2", "rfp_mat_rate3", "ntp_deg_rate", "nmp_deg_rate", "Km_ntp_rnap", "Km_ntp_ribo", "Km_aa_ribo", "kon_prom", "koff_prom_gfp", "koff_prom_rfp", "kon_mrna_gfp", "kon_mrna_rfp", "koff_mrna_gfp", "koff_mrna_rfp", "mrna_gfp_len", "mrna_rfp_len", "prot_gfp_len", "prot_rfp_len", "regen_rate", "Km_IInd_nrg_src", "Km_nmp_regen", "xylr_conc", "xylose_conc", "Kd", "Kx", "n", "m", "ribo_deg_rate",  "" };
const char* x_names[] = {"NMP", "NTP", "NTP_tx_gfp_track", "NTP_tx_rfp_track", "NTP_tl_gfp_track", "Ribo", "RNAP", "mRNA_RFP", "Ribo_mRNA_GFP", "Ribo_mRNA_RFP_elong", "Ribo_mRNA_GFP_elong", "RFP_i2", "RNAP_prom_GFP", "NTP_tl_rfp_track", "RFP_i1", "RNAP_prom_RFP", "GFP", "Ribo_mRNA_RFP", "IInd_nrg_src", "RFP", "RNAP_prom_RFP_elong", "mRNA_GFP",  "" };
const char* y_names[] = {"NTP_waste_track", "prom_GFP", "Mat_GFP", "prom_RFP", "RNAP_prom_GFP_elong", "AA", "Mat_RFP", "mRNA_GFP_sum", "mRNA_RFP_sum", "Ribo_per_mRNA_GFP_sum", "Ribo_per_mRNA_RFP_sum", "NTP_tx_total_track", "NTP_tl_total_track", "AA_tld_per_codon_txd", "RNAP_per_prom_GFP", "RNAP_per_prom_RFP",  "" };
const char* xc_names[] = {"NMP", "NTP", "NTP_tx_gfp_track", "NTP_tx_rfp_track", "NTP_tl_gfp_track", "Ribo", "RNAP", "mRNA_RFP", "Ribo_mRNA_GFP", "Ribo_mRNA_RFP_elong", "Ribo_mRNA_GFP_elong", "RFP_i2", "RNAP_prom_GFP", "NTP_tl_rfp_track", "RFP_i1", "RNAP_prom_RFP", "GFP", "Ribo_mRNA_RFP", "IInd_nrg_src", "RFP", "RNAP_prom_RFP_elong", "mRNA_GFP",  "" };
const char* pc_names[] = { "" };
const char* yc_names[] = {"NTP_waste_track", "prom_GFP", "Mat_GFP", "prom_RFP", "RNAP_prom_GFP_elong", "AA", "Mat_RFP", "mRNA_GFP_sum", "mRNA_RFP_sum", "Ribo_per_mRNA_GFP_sum", "Ribo_per_mRNA_RFP_sum", "NTP_tx_total_track", "NTP_tl_total_track", "AA_tld_per_codon_txd", "RNAP_per_prom_GFP", "RNAP_per_prom_RFP",  "" };
const char* dx_names[] = {"ODE NMP", "ODE NTP", "ODE NTP_tx_gfp_track", "ODE NTP_tx_rfp_track", "ODE NTP_tl_gfp_track", "ODE Ribo", "ODE RNAP", "ODE mRNA_RFP", "ODE Ribo_mRNA_GFP", "ODE Ribo_mRNA_RFP_elong", "ODE Ribo_mRNA_GFP_elong", "ODE RFP_i2", "ODE RNAP_prom_GFP", "ODE NTP_tl_rfp_track", "ODE RFP_i1", "ODE RNAP_prom_RFP", "ODE GFP", "ODE Ribo_mRNA_RFP", "ODE IInd_nrg_src", "ODE RFP", "ODE RNAP_prom_RFP_elong", "ODE mRNA_GFP",  "" };
const char* ct_names[] = {"CT NTP_waste_track", "CT prom_GFP", "CT Mat_GFP", "CT prom_RFP", "CT RNAP_prom_GFP_elong", "CT AA", "CT Mat_RFP", "CT mRNA_GFP_sum", "CT mRNA_RFP_sum", "CT Ribo_per_mRNA_GFP_sum", "CT Ribo_per_mRNA_RFP_sum", "CT NTP_tx_total_track", "CT NTP_tl_total_track", "CT AA_tld_per_codon_txd", "CT RNAP_per_prom_GFP", "CT RNAP_per_prom_RFP",  "" };
#endif // NAME_ARRAYS

#ifdef INITIAL
x[0] = 0;	//metabolite 'NMP': reactions
x[1] = 5000;	//metabolite 'NTP': reactions
x[2] = 0;	//metabolite 'NTP_tx_gfp_track': reactions
x[3] = 0;	//metabolite 'NTP_tx_rfp_track': reactions
x[4] = 0;	//metabolite 'NTP_tl_gfp_track': reactions
x[5] = 10;	//metabolite 'Ribo': reactions
x[6] = 1;	//metabolite 'RNAP': reactions
x[7] = 0;	//metabolite 'mRNA_RFP': reactions
x[8] = 0;	//metabolite 'Ribo_mRNA_GFP': reactions
x[9] = 0;	//metabolite 'Ribo_mRNA_RFP_elong': reactions
x[10] = 0;	//metabolite 'Ribo_mRNA_GFP_elong': reactions
x[11] = 0;	//metabolite 'RFP_i2': reactions
x[12] = 0;	//metabolite 'RNAP_prom_GFP': reactions
x[13] = 0;	//metabolite 'NTP_tl_rfp_track': reactions
x[14] = 0;	//metabolite 'RFP_i1': reactions
x[15] = 0;	//metabolite 'RNAP_prom_RFP': reactions
x[16] = 0;	//metabolite 'GFP': reactions
x[17] = 0;	//metabolite 'Ribo_mRNA_RFP': reactions
x[18] = 5000;	//metabolite 'IInd_nrg_src': reactions
x[19] = 0;	//metabolite 'RFP': reactions
x[20] = 0;	//metabolite 'RNAP_prom_RFP_elong': reactions
x[21] = 0;	//metabolite 'mRNA_GFP': reactions
#endif /* INITIAL */

#ifdef FIXED
ct[0] = 9999.999999999996;	//ct[0] conserved total for 'NTP_waste_track'
ct[1] = 1;	//ct[1] conserved total for 'prom_GFP'
ct[2] = 0;	//ct[2] conserved total for 'Mat_GFP'
ct[3] = 1;	//ct[3] conserved total for 'prom_RFP'
ct[4] = 1;	//ct[4] conserved total for 'RNAP_prom_GFP_elong'
ct[5] = 30000;	//ct[5] conserved total for 'AA'
ct[6] = 0;	//ct[6] conserved total for 'Mat_RFP'
ct[7] = 0; // manual add conserved total for mRNA_GFP_sum
ct[8] = 0; // manual add conserved total for mRNA_RFP_sum
ct[9] = 0; // manual add conserved total for Ribo_per_mRNA_GFP_sum
ct[10] = 0; // manual add conserved total for Ribo_per_mRNA_RFP_sum
ct[11] = 0; // manual add conserved total
ct[12] = 0; // manual add conserved total
ct[13] = 0; // manual add conserved total
ct[14] = 0;
ct[15] = 0;
p[0] = 1;	//compartment 'TestTube':fixed
p[1] = 10000;	//global quantity 'tl_rate_gfp':fixed
p[2] = 10000;	//global quantity 'tl_rate_rfp':fixed
p[3] = 100;	//global quantity 'tl_escape_rate_gfp':fixed
p[4] = 100;	//global quantity 'tl_escape_rate_rfp':fixed
p[5] = 1000;	//global quantity 'tx_rate':fixed
p[6] = 100;	//global quantity 'tx_escape_rate_gfp':fixed
p[7] = 100;	//global quantity 'tx_escape_rate_rfp':fixed
p[8] = 0.2;	//global quantity 'rna_deg_rate':fixed
p[9] = 0.06;	//global quantity 'gfp_mat_rate':fixed
p[10] = 0.01;	//global quantity 'rfp_mat_rate1':fixed
p[11] = 0.01;	//global quantity 'rfp_mat_rate2':fixed
p[12] = 0.01;	//global quantity 'rfp_mat_rate3':fixed
p[13] = 0.01;	//global quantity 'ntp_deg_rate':fixed
p[14] = 0.01;	//global quantity 'nmp_deg_rate':fixed
p[15] = 0.1;	//global quantity 'Km_ntp_rnap':fixed
p[16] = 10;	//global quantity 'Km_ntp_ribo':fixed
p[17] = 1000;	//global quantity 'Km_aa_ribo':fixed
p[18] = 1;	//global quantity 'kon_prom':fixed
p[19] = 4;	//global quantity 'koff_prom_gfp':fixed
p[20] = 4;	//global quantity 'koff_prom_rfp':fixed
p[21] = 400;	//global quantity 'kon_mrna_gfp':fixed
p[22] = 400;	//global quantity 'kon_mrna_rfp':fixed
p[23] = 100;	//global quantity 'koff_mrna_gfp':fixed
p[24] = 100;	//global quantity 'koff_mrna_rfp':fixed
p[25] = 1150;	//global quantity 'mrna_gfp_len':fixed
p[26] = 1102;	//global quantity 'mrna_rfp_len':fixed
p[27] = 245;	//global quantity 'prot_gfp_len':fixed
p[28] = 236;	//global quantity 'prot_rfp_len':fixed
p[29] = 60;	//global quantity 'regen_rate':fixed
p[30] = 800;	//global quantity 'Km_IInd_nrg_src':fixed
p[31] = 15;	//global quantity 'Km_nmp_regen':fixed
p[32] = 0;	//global quantity 'xylr_conc':fixed
p[33] = 0;	//global quantity 'xylose_conc':fixed
p[34] = 1;	//global quantity 'Kd':fixed
p[35] = 1;	//global quantity 'Kx':fixed
p[36] = 1;	//global quantity 'n':fixed
p[37] = 1;	//global quantity 'm':fixed
p[38] = 0.01;	//global quantity 'ribo_deg_rate':fixed
#endif /* FIXED */

#ifdef ASSIGNMENT
y[0] = ct[0]-1*x[1]-x[2]-1*x[3]-x[4]+9.73093e-14*x[7]+9.51977e-14*x[8]-1.03183e-13*x[11]-1*x[13]+7.85282e-14*x[14]+3.87607e-14*x[16]+9.79769e-14*x[17]-x[18]-3.22525e-13*x[19]+9.70836e-14*x[21];	//metabolite 'NTP_waste_track': reactions
y[1] = ct[1]-1.56418e-13*x[7]-1.63302e-13*x[8]-1*x[12]-1.56393e-13*x[17]-1.63326e-13*x[21];	//metabolite 'prom_GFP': reactions
y[2] = ct[2]+0.00204082*x[4]-5.23304e-14*x[7]-5.46764e-14*x[8]-x[16]-5.23198e-14*x[17]-5.46742e-14*x[21];	//metabolite 'Mat_GFP': reactions
y[3] = ct[3]+4.52636e-14*x[7]+4.72172e-14*x[8]-x[15]+4.52697e-14*x[17]+4.72113e-14*x[21];	//metabolite 'prom_RFP': reactions
y[4] = ct[4]-x[6]-1.96313e-13*x[7]-2.04949e-13*x[8]-1*x[12]-x[15]-1.96283e-13*x[17]-1*x[20]-2.04976e-13*x[21];	//metabolite 'RNAP_prom_GFP_elong': reactions
y[5] = ct[5]-0.5*x[4]+9.2662e-14*x[7]+9.04334e-14*x[8]-0.5*x[13]-2.99352e-14*x[14]+9.30142e-14*x[17]-3.6958e-14*x[19]+9.13116e-14*x[21];	//metabolite 'AA': reactions
y[6] = ct[6]-x[11]+0.00211864*x[13]-x[14]-x[19];	//metabolite 'Mat_RFP': reactions
y[7] = ct[7] + x[21] + x[8]; // manual add - sum of mRNA_GFP and Ribo_mRNA_GFP
y[8] = ct[8] + x[7] + x[17]; // manual add - sum of mRNA_RFP and Ribo_mRNA_RFP
y[9] =  y[7] < 1e-9 ? 0 : ( (x[8]  + x[10]) / y[7] );			// manual add - (Ribo_mRNA_GFP + Ribo_mRNA_GFP_elong) / mRNA_GFP_sum 
y[10] = y[8] < 1e-9 ? 0 : ( (x[17] + x[9] ) / y[8] );			// manual add - (Ribo_mRNA_RFP + Ribo_mRNA_RFP_elong) / mRNA_RFP_sum
y[11] = x[2] + x[3];	// manual add - NTP_tx_gfp_track + NTP_tx_rfp_track
y[12] = x[4] + x[13];	// manual add - NTP_tl_gfp_track + NTP_tl_rfp_track
y[13] = y[11] < 1e-5 ? 0 : (y[12]/2.0)/(y[11]/3.0);   // manual add - (NTP_tl_total_track/2.0) / (NTP_tx_total_track/3.0)
y[14] = (x[12] + y[1]) < 1e-9 ? 0 : (x[12] + y[4]) / (x[12] + y[1]);		//(x[12] + y[1]);//(x[12] + y[4]) / (x[12] + y[1]);	// manual add (RNAP_per_prom_GFP) - (RNAP_prom_GFP + RNAP_prom_GFP_elong) / (RNAP_prom_GFP + prom_GFP)
y[15] = (x[15] + y[3]) < 1e-9 ? 0 : (x[15] + x[20])/ (x[15] + y[3]);	//(x[15] + y[3]);//(x[15] + x[20]) / (x[15] + y[3]);	// manual add (RNAP_per_prom_RFP) - (RNAP_prom_RFP + RNAP_prom_RFP_elong) / (RNAP_prom_RFP + prom_RFP)
x_c[0] = x[0]/p[0];	//concentration of metabolite 'NMP': reactions
x_c[1] = x[1]/p[0];	//concentration of metabolite 'NTP': reactions
x_c[2] = x[2]/p[0];	//concentration of metabolite 'NTP_tx_gfp_track': reactions
x_c[3] = x[3]/p[0];	//concentration of metabolite 'NTP_tx_rfp_track': reactions
x_c[4] = x[4]/p[0];	//concentration of metabolite 'NTP_tl_gfp_track': reactions
x_c[5] = x[5]/p[0];	//concentration of metabolite 'Ribo': reactions
x_c[6] = x[6]/p[0];	//concentration of metabolite 'RNAP': reactions
x_c[7] = x[7]/p[0];	//concentration of metabolite 'mRNA_RFP': reactions
x_c[8] = x[8]/p[0];	//concentration of metabolite 'Ribo_mRNA_GFP': reactions
x_c[9] = x[9]/p[0];	//concentration of metabolite 'Ribo_mRNA_RFP_elong': reactions
x_c[10] = x[10]/p[0];	//concentration of metabolite 'Ribo_mRNA_GFP_elong': reactions
x_c[11] = x[11]/p[0];	//concentration of metabolite 'RFP_i2': reactions
x_c[12] = x[12]/p[0];	//concentration of metabolite 'RNAP_prom_GFP': reactions
x_c[13] = x[13]/p[0];	//concentration of metabolite 'NTP_tl_rfp_track': reactions
x_c[14] = x[14]/p[0];	//concentration of metabolite 'RFP_i1': reactions
x_c[15] = x[15]/p[0];	//concentration of metabolite 'RNAP_prom_RFP': reactions
x_c[16] = x[16]/p[0];	//concentration of metabolite 'GFP': reactions
x_c[17] = x[17]/p[0];	//concentration of metabolite 'Ribo_mRNA_RFP': reactions
x_c[18] = x[18]/p[0];	//concentration of metabolite 'IInd_nrg_src': reactions
x_c[19] = x[19]/p[0];	//concentration of metabolite 'RFP': reactions
x_c[20] = x[20]/p[0];	//concentration of metabolite 'RNAP_prom_RFP_elong': reactions
x_c[21] = x[21]/p[0];	//concentration of metabolite 'mRNA_GFP': reactions
y_c[0] = y[0]/p[0];	//concentration of metabolite 'NTP_waste_track': reactions
y_c[1] = y[1]/p[0];	//concentration of metabolite 'prom_GFP': reactions
y_c[2] = y[2]/p[0];	//concentration of metabolite 'Mat_GFP': reactions
y_c[3] = y[3]/p[0];	//concentration of metabolite 'prom_RFP': reactions
y_c[4] = y[4]/p[0];	//concentration of metabolite 'RNAP_prom_GFP_elong': reactions
y_c[5] = y[5]/p[0];	//concentration of metabolite 'AA': reactions
y_c[6] = y[6]/p[0];	//concentration of metabolite 'Mat_RFP': reactions
y_c[7] = y[7]/p[0];     // manual edit
y_c[8] = y[8]/p[0];     // manual edit
y_c[9] = y[9]/p[0];     // manual edit
y_c[10] = y[10]/p[0];     // manual edit
y_c[11] = y[11]/p[0];     // manual edit
y_c[12] = y[12]/p[0];     // manual edit
y_c[13] = y[13]/p[0];     // manual edit
y_c[14] = y[14]/p[0];     // manual edit
y_c[15] = y[15]/p[0];     // manual edit
#endif /* ASSIGNMENT */

#ifdef FUNCTIONS_HEADERS
double FunctionForMrna_gfp_deg_rxn(double sub_0, double param_0); 
double FunctionForMrna_rfp_deg_rxn(double sub_0, double param_0); 
double FunctionForTx_gfp_rxn(double param_0, double sub_0, double sub_1, double param_1, double param_2); 
double FunctionForTx_rfp_rxn(double param_0, double sub_0, double sub_1, double param_1, double param_2); 
double FunctionForTl_gfp_rxn(double sub_0, double param_0, double param_1, double sub_1, double sub_2, double param_2, double param_3); 
double FunctionForTl_rfp_rxn(double sub_0, double param_0, double param_1, double sub_1, double sub_2, double param_2, double param_3); 
double FunctionForRegen_rxn(double sub_0, double param_0, double param_1, double sub_1, double param_2); 
double FunctionForTx_gfp_escape_rxn(double param_0, double modif_0, double sub_0, double param_1); 
double FunctionForTx_rfp_escape_rxn(double param_0, double modif_0, double sub_0, double param_1); 
double FunctionForTl_gfp_escape_rxn(double param_0, double modif_0, double sub_0, double param_1); 
double FunctionForTl_rfp_escape_rxn(double param_0, double modif_0, double sub_0, double param_1); 
double FunctionForRibo_mrna_gfp_deg_rxn(double sub_0, double param_0); 
double FunctionForRibo_mrna_rfp_deg_rxn(double sub_0, double param_0); 
#endif /* FUNCTIONS_HEADERS */

#ifdef FUNCTIONS
double FunctionForMrna_gfp_deg_rxn(double sub_0, double param_0) 	//Function for mrna_gfp_deg_rxn
{return  param_0*sub_0;} 
double FunctionForMrna_rfp_deg_rxn(double sub_0, double param_0) 	//Function for mrna_rfp_deg_rxn
{return  param_0*sub_0;} 
double FunctionForTx_gfp_rxn(double param_0, double sub_0, double sub_1, double param_1, double param_2) 	//Function for tx_gfp_rxn
{return  param_2*sub_1*(sub_0/param_0/(param_1*(1.00000000000000000+sub_0/param_0)));} 
double FunctionForTx_rfp_rxn(double param_0, double sub_0, double sub_1, double param_1, double param_2) 	//Function for tx_rfp_rxn
{return  param_2*sub_1*(sub_0/param_0/(param_1*(1.00000000000000000+sub_0/param_0)));} 
double FunctionForTl_gfp_rxn(double sub_0, double param_0, double param_1, double sub_1, double sub_2, double param_2, double param_3) 	//Function for tl_gfp_rxn
{return  param_3*(sub_0/param_0*(sub_1/param_1)/(param_2*(1.00000000000000000+sub_0/param_0)*(1.00000000000000000+sub_1/param_1)))*sub_2;} 
double FunctionForTl_rfp_rxn(double sub_0, double param_0, double param_1, double sub_1, double sub_2, double param_2, double param_3) 	//Function for tl_rfp_rxn
{return  param_3*(sub_0/param_0*(sub_1/param_1)/(param_2*(1.00000000000000000+sub_0/param_0)*(1.00000000000000000+sub_1/param_1)))*sub_2;} 
double FunctionForRegen_rxn(double sub_0, double param_0, double param_1, double sub_1, double param_2) 	//Function for regen_rxn
{return  param_2*(sub_1/param_1*(sub_0/param_0)/((1.00000000000000000+sub_1/param_1)*(1.00000000000000000+sub_0/param_0)));} 
double FunctionForTx_gfp_escape_rxn(double param_0, double modif_0, double sub_0, double param_1) 	//Function for tx_gfp_escape_rxn
{return  param_1*sub_0*(modif_0/param_0/(1.00000000000000000+modif_0/param_0));} 
double FunctionForTx_rfp_escape_rxn(double param_0, double modif_0, double sub_0, double param_1) 	//Function for tx_rfp_escape_rxn
{return  param_1*sub_0*(modif_0/param_0/(1.00000000000000000+modif_0/param_0));} 
double FunctionForTl_gfp_escape_rxn(double param_0, double modif_0, double sub_0, double param_1) 	//Function for tl_gfp_escape_rxn
{return  param_1*sub_0*(modif_0/param_0/(1.00000000000000000+modif_0/param_0));} 
double FunctionForTl_rfp_escape_rxn(double param_0, double modif_0, double sub_0, double param_1) 	//Function for tl_rfp_escape_rxn
{return  param_1*sub_0*(modif_0/param_0/(1.00000000000000000+modif_0/param_0));} 
double FunctionForRibo_mrna_gfp_deg_rxn(double sub_0, double param_0) 	//Function for ribo_mrna_gfp_deg_rxn
{return  param_0*sub_0;} 
double FunctionForRibo_mrna_rfp_deg_rxn(double sub_0, double param_0) 	//Function for ribo_mrna_rfp_deg_rxn
{return  param_0*sub_0;} 
#endif /* FUNCTIONS */

#ifdef ODEs
dx[0] = 1150*FunctionForMrna_gfp_deg_rxn(x_c[21], p[8])*p[0]+1102*FunctionForMrna_rfp_deg_rxn(x_c[7], p[8])*p[0]+(p[13] * x_c[1]) *p[0]-(p[14] * x_c[0]) *p[0]+490*FunctionForTl_gfp_rxn(y_c[5], p[17], p[16], x_c[1], x_c[10], p[27], p[1])*p[0]+472*FunctionForTl_rfp_rxn(y_c[5], p[17], p[16], x_c[1], x_c[9], p[28], p[2])*p[0]-FunctionForRegen_rxn(x_c[18], p[30], p[31], x_c[0], p[29])*p[0]+1150*FunctionForRibo_mrna_gfp_deg_rxn(x_c[8], p[8])*p[0]+1102*FunctionForRibo_mrna_rfp_deg_rxn(x_c[17], p[8])*p[0];	// 
dx[1] = -(p[13] * x_c[1]) *p[0]-1150*FunctionForTx_gfp_rxn(p[15], x_c[1], y_c[4], p[25], p[5])*p[0]-1102*FunctionForTx_rfp_rxn(p[15], x_c[1], x_c[20], p[26], p[5])*p[0]-490*FunctionForTl_gfp_rxn(y_c[5], p[17], p[16], x_c[1], x_c[10], p[27], p[1])*p[0]-472*FunctionForTl_rfp_rxn(y_c[5], p[17], p[16], x_c[1], x_c[9], p[28], p[2])*p[0]+FunctionForRegen_rxn(x_c[18], p[30], p[31], x_c[0], p[29])*p[0];	// 
dx[2] = 1150*FunctionForTx_gfp_rxn(p[15], x_c[1], y_c[4], p[25], p[5])*p[0];	// 
dx[3] = 1102*FunctionForTx_rfp_rxn(p[15], x_c[1], x_c[20], p[26], p[5])*p[0];	// 
dx[4] = 490*FunctionForTl_gfp_rxn(y_c[5], p[17], p[16], x_c[1], x_c[10], p[27], p[1])*p[0];	// 
dx[5] = -(p[21] * x_c[21] * x_c[5]) *p[0]-(p[22] * x_c[7] * x_c[5]) *p[0]+(p[23] * x_c[8]) *p[0]+(p[24] * x_c[17]) *p[0]+FunctionForTl_gfp_rxn(y_c[5], p[17], p[16], x_c[1], x_c[10], p[27], p[1])*p[0]+FunctionForTl_rfp_rxn(y_c[5], p[17], p[16], x_c[1], x_c[9], p[28], p[2])*p[0]+FunctionForRibo_mrna_gfp_deg_rxn(x_c[8], p[8])*p[0]+FunctionForRibo_mrna_rfp_deg_rxn(x_c[17], p[8])*p[0]+(p[8] * x_c[10]) *p[0]+(p[8] * x_c[9]) *p[0]-(p[38] * x_c[5]) *p[0];	// 
dx[6] = -(p[18] * y_c[1] * x_c[6]) *p[0]-(p[18] * y_c[3] * x_c[6]) *p[0]+(p[19] * x_c[12]) *p[0]+(p[20] * x_c[15]) *p[0]+FunctionForTx_gfp_rxn(p[15], x_c[1], y_c[4], p[25], p[5])*p[0]+FunctionForTx_rfp_rxn(p[15], x_c[1], x_c[20], p[26], p[5])*p[0];	// 
dx[7] = -FunctionForMrna_rfp_deg_rxn(x_c[7], p[8])*p[0]-(p[22] * x_c[7] * x_c[5]) *p[0]+(p[24] * x_c[17]) *p[0]+FunctionForTx_rfp_rxn(p[15], x_c[1], x_c[20], p[26], p[5])*p[0]+FunctionForTl_rfp_escape_rxn(p[16], x_c[1], x_c[17], p[4])*p[0]+(p[38] * x_c[17]) *p[0];	// 
dx[8] = (p[21] * x_c[21] * x_c[5]) *p[0]-(p[23] * x_c[8]) *p[0]-FunctionForTl_gfp_escape_rxn(p[16], x_c[1], x_c[8], p[3])*p[0]-FunctionForRibo_mrna_gfp_deg_rxn(x_c[8], p[8])*p[0]-(p[38] * x_c[8]) *p[0];	// 
dx[9] = -FunctionForTl_rfp_rxn(y_c[5], p[17], p[16], x_c[1], x_c[9], p[28], p[2])*p[0]+FunctionForTl_rfp_escape_rxn(p[16], x_c[1], x_c[17], p[4])*p[0]-(p[8] * x_c[9]) *p[0]-(p[38] * x_c[9]) *p[0];	// 
dx[10] = -FunctionForTl_gfp_rxn(y_c[5], p[17], p[16], x_c[1], x_c[10], p[27], p[1])*p[0]+FunctionForTl_gfp_escape_rxn(p[16], x_c[1], x_c[8], p[3])*p[0]-(p[8] * x_c[10]) *p[0]-(p[38] * x_c[10]) *p[0];	// 
dx[11] = (p[11] * x_c[14]) *p[0]-(p[12] * x_c[11]) *p[0];	// 
dx[12] = (p[18] * y_c[1] * x_c[6]) *p[0]-(p[19] * x_c[12]) *p[0]-FunctionForTx_gfp_escape_rxn(p[15], x_c[1], x_c[12], p[6])*p[0];	// 
dx[13] = 472*FunctionForTl_rfp_rxn(y_c[5], p[17], p[16], x_c[1], x_c[9], p[28], p[2])*p[0];	// 
dx[14] = (p[10] * x_c[19]) *p[0]-(p[11] * x_c[14]) *p[0];	// 
dx[15] = (p[18] * y_c[3] * x_c[6]) *p[0]-(p[20] * x_c[15]) *p[0]-FunctionForTx_rfp_escape_rxn(p[15], x_c[1], x_c[15], p[7])*p[0];	// 
dx[16] = -(p[9] * x_c[16]) *p[0]+FunctionForTl_gfp_rxn(y_c[5], p[17], p[16], x_c[1], x_c[10], p[27], p[1])*p[0];	// 
dx[17] = (p[22] * x_c[7] * x_c[5]) *p[0]-(p[24] * x_c[17]) *p[0]-FunctionForTl_rfp_escape_rxn(p[16], x_c[1], x_c[17], p[4])*p[0]-FunctionForRibo_mrna_rfp_deg_rxn(x_c[17], p[8])*p[0]-(p[38] * x_c[17]) *p[0];	// 
dx[18] = -FunctionForRegen_rxn(x_c[18], p[30], p[31], x_c[0], p[29])*p[0];	// 
dx[19] = -(p[10] * x_c[19]) *p[0]+FunctionForTl_rfp_rxn(y_c[5], p[17], p[16], x_c[1], x_c[9], p[28], p[2])*p[0];	// 
dx[20] = -FunctionForTx_rfp_rxn(p[15], x_c[1], x_c[20], p[26], p[5])*p[0]+FunctionForTx_rfp_escape_rxn(p[15], x_c[1], x_c[15], p[7])*p[0];	// 
dx[21] = -FunctionForMrna_gfp_deg_rxn(x_c[21], p[8])*p[0]-(p[21] * x_c[21] * x_c[5]) *p[0]+(p[23] * x_c[8]) *p[0]+FunctionForTx_gfp_rxn(p[15], x_c[1], y_c[4], p[25], p[5])*p[0]+FunctionForTl_gfp_escape_rxn(p[16], x_c[1], x_c[8], p[3])*p[0]+(p[38] * x_c[8]) *p[0];	// 
#endif /* ODEs */

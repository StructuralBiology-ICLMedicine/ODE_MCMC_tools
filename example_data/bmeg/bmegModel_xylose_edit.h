#ifdef SIZE_DEFINITIONS
#define N_METABS 14
#define N_ODE_METABS 0
#define N_INDEP_METABS 11
#define N_COMPARTMENTS 1
#define N_GLOBAL_PARAMS 27
#define N_KIN_PARAMS 0
#define N_REACTIONS 18

#define N_ARRAY_SIZE_P  28	// number of parameters
#define N_ARRAY_SIZE_X  11	// number of initials
#define N_ARRAY_SIZE_Y  4	//	3	// number of assigned elements
#define N_ARRAY_SIZE_XC 11	// number of x concentration
#define N_ARRAY_SIZE_PC 0	// number of p concentration
#define N_ARRAY_SIZE_YC 4	//	3	// number of y concentration
#define N_ARRAY_SIZE_DX 11	// number of ODEs 
#define N_ARRAY_SIZE_CT 4	//	3	// number of conserved totals

#endif // SIZE_DEFINITIONS

#ifdef TIME
#define T  <set here a user name for the time variable> 
#endif // TIME

#ifdef NAME_ARRAYS
const char* p_names[] = {"TestTube", "tl_rate_gfp", "tl_escape_rate_gfp", "tx_rate", "tx_escape_rate_gfp", "rna_deg_rate_gfp", "gfp_mat_rate", "ntp_deg_rate", "nmp_deg_rate", "Km_ntp_rnap", "Km_ntp_ribo", "Km_aa_ribo", "kon_prom", "koff_prom_gfp", "kon_mrna_gfp", "koff_mrna_gfp", "mrna_gfp_len", "prot_gfp_len", "regen_rate", "Km_IInd_nrg_src", "Km_nmp_regen", "xylr_conc", "xylose_conc", "Kd", "Kx", "n", "m", "ribo_deg_rate",  "" };
const char* x_names[] = {"NMP", "NTP", "AA", "mRNA_GFP", "Ribo", "prom_GFP", "Ribo_mRNA_GFP_elong", "GFP", "IInd_nrg_src", "RNAP", "Ribo_mRNA_GFP",  "" };
const char* y_names[] = {"Mat_GFP", "RNAP_prom_GFP_elong", "RNAP_prom_GFP", "mRNA_GFP_sum",  "" };
const char* xc_names[] = {"NMP", "NTP", "AA", "mRNA_GFP", "Ribo", "prom_GFP", "Ribo_mRNA_GFP_elong", "GFP", "IInd_nrg_src", "RNAP", "Ribo_mRNA_GFP",  "" };
const char* pc_names[] = { "" };
const char* yc_names[] = {"Mat_GFP", "RNAP_prom_GFP_elong", "RNAP_prom_GFP", "mRNA_GFP_sum",  "" };
const char* dx_names[] = {"ODE NMP", "ODE NTP", "ODE AA", "ODE mRNA_GFP", "ODE Ribo", "ODE prom_GFP", "ODE Ribo_mRNA_GFP_elong", "ODE GFP", "ODE IInd_nrg_src", "ODE RNAP", "ODE Ribo_mRNA_GFP",  "" };
const char* ct_names[] = {"CT Mat_GFP", "CT RNAP_prom_GFP_elong", "CT RNAP_prom_GFP", "CT mRNA_GFP_sum", "" };
#endif // NAME_ARRAYS

#ifdef INITIAL
x[0] = 0;	//metabolite 'NMP': reactions
x[1] = 5000;	//metabolite 'NTP': reactions
x[2] = 30000;	//metabolite 'AA': reactions
x[3] = 0;	//metabolite 'mRNA_GFP': reactions
x[4] = 10;	//metabolite 'Ribo': reactions
x[5] = 1;	//metabolite 'prom_GFP': reactions
x[6] = 0;	//metabolite 'Ribo_mRNA_GFP_elong': reactions
x[7] = 0;	//metabolite 'GFP': reactions
x[8] = 5000;	//metabolite 'IInd_nrg_src': reactions
x[9] = 1;	//metabolite 'RNAP': reactions
x[10] = 0;	//metabolite 'Ribo_mRNA_GFP': reactions
#endif /* INITIAL */

#ifdef FIXED
ct[0] = 122.4489795918366;	//ct[0] conserved total for 'Mat_GFP'
ct[1] = -1.110223024625157e-16;	//ct[1] conserved total for 'RNAP_prom_GFP_elong'
ct[2] = 1;	//ct[2] conserved total for 'RNAP_prom_GFP'
ct[3] = 0; // manual add conserved total for mRNA_GFP_sum
p[0] = 1;	//compartment 'TestTube':fixed
p[1] = 10000;	//global quantity 'tl_rate_gfp':fixed
p[2] = 100;	//global quantity 'tl_escape_rate_gfp':fixed
p[3] = 1000;	//global quantity 'tx_rate':fixed
p[4] = 100;	//global quantity 'tx_escape_rate_gfp':fixed
p[5] = 0.2;	//global quantity 'rna_deg_rate_gfp':fixed
p[6] = 0.06;	//global quantity 'gfp_mat_rate':fixed
p[7] = 0.01;	//global quantity 'ntp_deg_rate':fixed
p[8] = 0.01;	//global quantity 'nmp_deg_rate':fixed
p[9] = 0.1;	//global quantity 'Km_ntp_rnap':fixed
p[10] = 10;	//global quantity 'Km_ntp_ribo':fixed
p[11] = 1000;	//global quantity 'Km_aa_ribo':fixed
p[12] = 1;	//global quantity 'kon_prom':fixed
p[13] = 4;	//global quantity 'koff_prom_gfp':fixed
p[14] = 400;	//global quantity 'kon_mrna_gfp':fixed
p[15] = 100;	//global quantity 'koff_mrna_gfp':fixed
p[16] = 1150;	//global quantity 'mrna_gfp_len':fixed
p[17] = 245;	//global quantity 'prot_gfp_len':fixed
p[18] = 60;	//global quantity 'regen_rate':fixed
p[19] = 800;	//global quantity 'Km_IInd_nrg_src':fixed
p[20] = 15;	//global quantity 'Km_nmp_regen':fixed
p[21] = 0;	//global quantity 'xylr_conc':fixed
p[22] = 0;	//global quantity 'xylose_conc':fixed
p[23] = 1;	//global quantity 'Kd':fixed
p[24] = 1;	//global quantity 'Kx':fixed
p[25] = 1;	//global quantity 'n':fixed
p[26] = 1;	//global quantity 'm':fixed
p[27] = 0.01;	//global quantity 'ribo_deg_rate':fixed
#endif /* FIXED */

#ifdef ASSIGNMENT
y[0] = ct[0]-0.00408163*x[2]+3.16218e-13*x[3]-1*x[7]+3.16223e-13*x[10];	//metabolite 'Mat_GFP': reactions
y[1] = ct[1]-7.79484e-14*x[3]+x[5]-x[9]-7.79454e-14*x[10];	//metabolite 'RNAP_prom_GFP_elong': reactions
y[2] = ct[2]-x[5];	//metabolite 'RNAP_prom_GFP': reactionsi
y[3] = ct[3] + x[3] + x[10]; // manual add - sum of mRNA_GFP and Ribo_mRNA_GFP
x_c[0] = x[0]/p[0];	//concentration of metabolite 'NMP': reactions
x_c[1] = x[1]/p[0];	//concentration of metabolite 'NTP': reactions
x_c[2] = x[2]/p[0];	//concentration of metabolite 'AA': reactions
x_c[3] = x[3]/p[0];	//concentration of metabolite 'mRNA_GFP': reactions
x_c[4] = x[4]/p[0];	//concentration of metabolite 'Ribo': reactions
x_c[5] = x[5]/p[0];	//concentration of metabolite 'prom_GFP': reactions
x_c[6] = x[6]/p[0];	//concentration of metabolite 'Ribo_mRNA_GFP_elong': reactions
x_c[7] = x[7]/p[0];	//concentration of metabolite 'GFP': reactions
x_c[8] = x[8]/p[0];	//concentration of metabolite 'IInd_nrg_src': reactions
x_c[9] = x[9]/p[0];	//concentration of metabolite 'RNAP': reactions
x_c[10] = x[10]/p[0];	//concentration of metabolite 'Ribo_mRNA_GFP': reactions
y_c[0] = y[0]/p[0];	//concentration of metabolite 'Mat_GFP': reactions
y_c[1] = y[1]/p[0];	//concentration of metabolite 'RNAP_prom_GFP_elong': reactions
y_c[2] = y[2]/p[0];	//concentration of metabolite 'RNAP_prom_GFP': reactions
y_c[3] = y[3]/p[0];	// manual edit
#endif /* ASSIGNMENT */

#ifdef FUNCTIONS_HEADERS
double FunctionForMrna_gfp_deg_rxn(double sub_0, double param_0); 
double FunctionForTx_gfp_rxn_3(double param_0, double sub_0, double sub_1, double param_1, double param_2); 
double FunctionForTl_gfp_rxn_2(double sub_0, double param_0, double param_1, double sub_1, double sub_2, double param_2, double param_3); 
double FunctionForRegen_rxn(double sub_0, double param_0, double param_1, double sub_1, double param_2); 
double FunctionForTx_gfp_escape_rxn(double param_0, double modif_0, double sub_0, double param_1); 
double FunctionForTl_gfp_escape_rxn(double param_0, double modif_0, double sub_0, double param_1); 
double FunctionForRibo_mrna_gfp_deg_rxn(double sub_0, double param_0); 
#endif /* FUNCTIONS_HEADERS */

#ifdef FUNCTIONS
double FunctionForMrna_gfp_deg_rxn(double sub_0, double param_0) 	//Function for mrna_gfp_deg_rxn
{return  param_0*sub_0;} 
double FunctionForTx_gfp_rxn_3(double param_0, double sub_0, double sub_1, double param_1, double param_2) 	//Function for tx_gfp_rxn_3
{return  param_2*sub_1*(sub_0/param_0/(param_1*(1.00000000000000000+sub_0/param_0)));} 
double FunctionForTl_gfp_rxn_2(double sub_0, double param_0, double param_1, double sub_1, double sub_2, double param_2, double param_3) 	//Function for tl_gfp_rxn_2
{return  param_3*(sub_0/param_0*(sub_1/param_1)/(param_2*(1.00000000000000000+sub_0/param_0)*(1.00000000000000000+sub_1/param_1)))*sub_2;} 
double FunctionForRegen_rxn(double sub_0, double param_0, double param_1, double sub_1, double param_2) 	//Function for regen_rxn
{return  param_2*(sub_1/param_1*(sub_0/param_0)/((1.00000000000000000+sub_1/param_1)*(1.00000000000000000+sub_0/param_0)));} 
double FunctionForTx_gfp_escape_rxn(double param_0, double modif_0, double sub_0, double param_1) 	//Function for tx_gfp_escape_rxn
{return  param_1*sub_0*(modif_0/param_0/(1.00000000000000000+modif_0/param_0));} 
double FunctionForTl_gfp_escape_rxn(double param_0, double modif_0, double sub_0, double param_1) 	//Function for tl_gfp_escape_rxn
{return  param_1*sub_0*(modif_0/param_0/(1.00000000000000000+modif_0/param_0));} 
double FunctionForRibo_mrna_gfp_deg_rxn(double sub_0, double param_0) 	//Function for ribo_mrna_gfp_deg_rxn
{return  param_0*sub_0;} 
#endif /* FUNCTIONS */

#ifdef ODEs
dx[0] = 1150*FunctionForMrna_gfp_deg_rxn(x_c[3], p[5])*p[0]+(p[7] * x_c[1]) *p[0]-(p[8] * x_c[0]) *p[0]+490*FunctionForTl_gfp_rxn_2(x_c[2], p[11], p[10], x_c[1], x_c[6], p[17], p[1])*p[0]-FunctionForRegen_rxn(x_c[8], p[19], p[20], x_c[0], p[18])*p[0]+1150*FunctionForRibo_mrna_gfp_deg_rxn(x_c[10], p[5])*p[0];	// 
dx[1] = -(p[7] * x_c[1]) *p[0]-1150*FunctionForTx_gfp_rxn_3(p[9], x_c[1], y_c[1], p[16], p[3])*p[0]-490*FunctionForTl_gfp_rxn_2(x_c[2], p[11], p[10], x_c[1], x_c[6], p[17], p[1])*p[0]+FunctionForRegen_rxn(x_c[8], p[19], p[20], x_c[0], p[18])*p[0];	// 
dx[2] = -245*FunctionForTl_gfp_rxn_2(x_c[2], p[11], p[10], x_c[1], x_c[6], p[17], p[1])*p[0];	// 
dx[3] = -FunctionForMrna_gfp_deg_rxn(x_c[3], p[5])*p[0]-(p[14] * x_c[3] * x_c[4]) *p[0]+(p[15] * x_c[10]) *p[0]+FunctionForTx_gfp_rxn_3(p[9], x_c[1], y_c[1], p[16], p[3])*p[0]+FunctionForTl_gfp_escape_rxn(p[10], x_c[1], x_c[10], p[2])*p[0]+(p[27] * x_c[10]) *p[0];	// 
dx[4] = -(p[14] * x_c[3] * x_c[4]) *p[0]+(p[15] * x_c[10]) *p[0]+FunctionForTl_gfp_rxn_2(x_c[2], p[11], p[10], x_c[1], x_c[6], p[17], p[1])*p[0]+FunctionForRibo_mrna_gfp_deg_rxn(x_c[10], p[5])*p[0]+(p[5] * x_c[6]) *p[0]-(p[27] * x_c[4]) *p[0];	// 
dx[5] = -(p[12] * x_c[5] * x_c[9]) *p[0]+(p[13] * y_c[2]) *p[0]+FunctionForTx_gfp_escape_rxn(p[9], x_c[1], y_c[2], p[4])*p[0];	// 
dx[6] = -FunctionForTl_gfp_rxn_2(x_c[2], p[11], p[10], x_c[1], x_c[6], p[17], p[1])*p[0]+FunctionForTl_gfp_escape_rxn(p[10], x_c[1], x_c[10], p[2])*p[0]-(p[5] * x_c[6]) *p[0]-(p[27] * x_c[6]) *p[0];	// 
dx[7] = -(p[6] * x_c[7]) *p[0]+FunctionForTl_gfp_rxn_2(x_c[2], p[11], p[10], x_c[1], x_c[6], p[17], p[1])*p[0];	// 
dx[8] = -FunctionForRegen_rxn(x_c[8], p[19], p[20], x_c[0], p[18])*p[0];	// 
dx[9] = -(p[12] * x_c[5] * x_c[9]) *p[0]+(p[13] * y_c[2]) *p[0]+FunctionForTx_gfp_rxn_3(p[9], x_c[1], y_c[1], p[16], p[3])*p[0];	// 
dx[10] = (p[14] * x_c[3] * x_c[4]) *p[0]-(p[15] * x_c[10]) *p[0]-FunctionForTl_gfp_escape_rxn(p[10], x_c[1], x_c[10], p[2])*p[0]-FunctionForRibo_mrna_gfp_deg_rxn(x_c[10], p[5])*p[0]-(p[27] * x_c[10]) *p[0];	// 
#endif /* ODEs */


#ifdef INITIAL_ASSIGN
const double Kx = p[24]; // edit
const double Kd = p[23]; // edit
const double m = p[26];  //edit
const double n = p[25];  //edit

//if (Kx > 0 && Kd > 0 && m > 0 && n > 0){
// initialise GFP conc using hill function

// must change in fact ct[3] and y[3] and y_c[3]
//y[3] = ct[3]-1*x[9];  //metabolite 'prom_GFP': reactions
const double init_prom_GFP = x_c[5];	//y_c[3];    //x_c[9]     ;//OLD //y_c[3];//x_c[7];
const double xylr_conc = p[21]; // edit
const double xylose_conc = p[22]; //edit

/*
std::cerr << "Kx: " << Kx << std::endl;
std::cerr << "Kd: " << Kd << std::endl;
std::cerr << "m: " << m << std::endl;
std::cerr << "n: " << n << std::endl;
*/

const double free_xylr = (xylr_conc)*(1/( std::pow(((xylose_conc)/Kx),m)  + 1));

/*
std::cerr << "INFO: free_xylr " << free_xylr << std::endl;
std::cerr << "INFO:  xylr_conc " << xylr_conc
                << " xylose_conc " << xylose_conc
                << std::endl;
 */

const double new_prom_GFP = (init_prom_GFP)* (1/( std::pow(((free_xylr)/Kd), n) + 1));


if (Kx > 0 && Kd > 0 && m > 0 && n > 0 && xylr_conc > 0  ){

        //std::cerr << "ERROR: need to fix this part\n" << std::endl;

        x[5] = new_prom_GFP*p[0];       //ct[3] = (new_prom_GFP*p[0]) + 1*x[9];
        //y[3] = ct[3]-1*x[9];
        x_c[5] = x[5] / p[0];
       //std::cerr << "INFO:\tchanging prom_GFP from " << init_prom_GFP << " to " << new_prom_GFP << " x_c[5] " << x_c[5] << std::endl;
}
else {
        //std::cerr << "INFO:\tNOT changing prom_GFP from " << init_prom_GFP << " to " << new_prom_GFP << std::endl;
}
//}
#endif /* INITIAL_ASSIGN */

//#ifndef BMEG_MODEL_H_
//#define BMEG_MODEL_H_

//#include "eigenmvn.h"
#include <string>
#include <vector>
#include <cmath>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include "common_model_defs.h"
#include "mcmc_interfaces.h"

namespace MODEL_CLASS_NAMESPACE {

#define SIZE_DEFINITIONS
#include MODEL_DEF_HEADER_FILE
#undef SIZE_DEFINITIONS


template<typename T, size_t N>
T * carr_end(T (&ra)[N])
{
    return ra + N;
}

class MODEL_CLASS_NAME : public MCMC::ode_model_system_interface {

private:


#define FUNCTIONS
#include MODEL_DEF_HEADER_FILE
#undef FUNCTIONS


public:
	
	
	virtual std::shared_ptr< MCMC::ode_model_system_interface > clone() override{
		std::shared_ptr< MCMC::ode_model_system_interface > new_ptr(new MODEL_CLASS_NAME(*this));
		return new_ptr;
	}


    /*
    // using boost:array
    typedef boost::array< double , N_ARRAY_SIZE_P > para_type;
    typedef boost::array< double , N_ARRAY_SIZE_X > state_type;
    typedef boost::array< double , N_ARRAY_SIZE_Y > state_y_type;
    typedef boost::array< double , N_ARRAY_SIZE_CT > CT_type;

    //! joint type for use in sampling both parameters and initial conditions
    typedef boost::array< double , N_ARRAY_SIZE_P+N_ARRAY_SIZE_X+N_ARRAY_SIZE_Y > joint_para_x_y_state_type;

    //! joint type for use in sampling both species types X and Y
    typedef boost::array< double , N_ARRAY_SIZE_X+N_ARRAY_SIZE_Y > joint_x_y_state_type;
    */


    // using boost::numeric::ublas::vector
    typedef boost::numeric::ublas::vector< double> para_type;
    static para_type new_para_type(){ return para_type(N_ARRAY_SIZE_P, 0); }
    typedef boost::numeric::ublas::vector< double> state_type;
    static state_type static_new_state_type() {
        return state_type(N_ARRAY_SIZE_X, 0);
    }
    virtual state_type new_state_type() const override{
        return static_new_state_type(); //state_type(N_ARRAY_SIZE_X, 0);
    }
    typedef boost::numeric::ublas::vector< double > state_y_type;
    static state_y_type new_state_y_type(){ return state_y_type(N_ARRAY_SIZE_Y, 0); }
    typedef boost::numeric::ublas::vector< double > CT_type;
    static CT_type new_CT_type(){ return CT_type(N_ARRAY_SIZE_CT, 0); }

    //! joint type for use in sampling both parameters and initial conditions
    typedef boost::numeric::ublas::vector< double   > joint_para_x_y_state_type;
    static joint_para_x_y_state_type new_joint_para_x_y_state_type(){ return joint_para_x_y_state_type(N_ARRAY_SIZE_P+N_ARRAY_SIZE_X+N_ARRAY_SIZE_Y, 0); }

    //! synonym for joint_para_x_y_state_type
    typedef joint_para_x_y_state_type joint_para_and_state_type;
    virtual joint_para_x_y_state_type new_joint_para_and_state_type() const override{
        return joint_para_x_y_state_type(N_ARRAY_SIZE_P+N_ARRAY_SIZE_X+N_ARRAY_SIZE_Y, 0);
    }

    virtual unsigned int get_dims() const override {
    	return N_ARRAY_SIZE_P+N_ARRAY_SIZE_X+N_ARRAY_SIZE_Y;
    }

    //! joint type for use in sampling both species types X and Y
    typedef boost::numeric::ublas::vector< double > joint_x_y_state_type;
    static joint_x_y_state_type new_joint_x_y_state_type(){ return joint_x_y_state_type(N_ARRAY_SIZE_X+N_ARRAY_SIZE_Y, 0); }

    typedef boost::numeric::ublas::matrix< double > state_matrix_type;
    virtual state_matrix_type new_state_matrix_type() const override{
        return state_matrix_type(N_ARRAY_SIZE_X, N_ARRAY_SIZE_X, 0);
    }


    para_type p;
    state_type x_c,
        x;
    state_y_type y_c,
        y;
    CT_type ct;

	const static std::vector<std::string> p_name_vec,
        x_name_vec,
        y_name_vec,
        xc_name_vec,
        pc_name_vec,
        yc_name_vec,
        dx_name_vec,
        ct_name_vec;


    static joint_x_y_state_type get_joint_x_y_state(const state_type x_, const state_y_type y_){
        joint_x_y_state_type joint = new_joint_x_y_state_type();
        std::copy(x_.begin(), x_.end(), joint.begin());
        std::copy(y_.begin(), y_.end(), joint.begin() + N_ARRAY_SIZE_X );
        return joint;
    }

    static state_type get_x_array_from_joint_x_y_state(joint_x_y_state_type joint){
        state_type x_ = static_new_state_type();
        std::copy(joint.begin(), joint.begin() + N_ARRAY_SIZE_X  , x_.begin());
        return x_;
    }

    static state_type get_y_array_from_joint_x_y_state(joint_x_y_state_type joint){
        state_type y_ = new_state_y_type();
        std::copy(joint.begin() + N_ARRAY_SIZE_X, joint.begin() + N_ARRAY_SIZE_X + N_ARRAY_SIZE_Y  , y_.begin());
        return y_;
    }


    static joint_para_x_y_state_type get_joint_para_x_y_state(const para_type p_, const state_type x_, const state_y_type y_){
        joint_para_x_y_state_type joint = new_joint_para_x_y_state_type();
        std::copy(p_.begin(), p_.end(), joint.begin());
        std::copy(x_.begin(), x_.end(), joint.begin() + N_ARRAY_SIZE_P );
        std::copy(y_.begin(), y_.end(), joint.begin() + N_ARRAY_SIZE_P + N_ARRAY_SIZE_X );
        return joint;
    }

    std::vector<std::string> get_joint_para_and_state_name_vec() const override{ // was called get_joint_para_x_y_state_name_vec()
        std::vector<std::string> joint_name_vec;
        joint_name_vec.reserve(N_ARRAY_SIZE_P + N_ARRAY_SIZE_X + N_ARRAY_SIZE_Y);
        joint_name_vec.insert(joint_name_vec.end(), p_name_vec.begin(), p_name_vec.begin()+N_ARRAY_SIZE_P);
        joint_name_vec.insert(joint_name_vec.end(), x_name_vec.begin(), x_name_vec.begin()+N_ARRAY_SIZE_X);
        joint_name_vec.insert(joint_name_vec.end(), y_name_vec.begin(), y_name_vec.begin()+N_ARRAY_SIZE_Y);

        return joint_name_vec;
    }

    static std::vector<std::string> get_joint_state_name_vec(){ // was called get_joint_para_x_y_state_name_vec()
            std::vector<std::string> joint_name_vec;
            joint_name_vec.reserve( N_ARRAY_SIZE_X + N_ARRAY_SIZE_Y);
            //joint_name_vec.insert(joint_name_vec.end(), p_name_vec.begin(), p_name_vec.begin()+N_ARRAY_SIZE_P);
            joint_name_vec.insert(joint_name_vec.end(), x_name_vec.begin(), x_name_vec.begin()+N_ARRAY_SIZE_X);
            joint_name_vec.insert(joint_name_vec.end(), y_name_vec.begin(), y_name_vec.begin()+N_ARRAY_SIZE_Y);

            return joint_name_vec;
        }

    static para_type get_p_array(joint_para_x_y_state_type joint){
        para_type p_ = new_para_type();
        std::copy(joint.begin(), joint.begin() + N_ARRAY_SIZE_P  , p_.begin());
        return p_;
    }


    static state_type get_x_array(joint_para_x_y_state_type joint){
        state_type x_ = static_new_state_type();
        std::copy(joint.begin() + N_ARRAY_SIZE_P, joint.begin() + N_ARRAY_SIZE_P + N_ARRAY_SIZE_X  , x_.begin());
        return x_;
    }

    static state_y_type get_y_array(joint_para_x_y_state_type joint){
        state_y_type y_ = new_state_y_type();
        std::copy(joint.begin() + N_ARRAY_SIZE_P + N_ARRAY_SIZE_X, joint.begin() + N_ARRAY_SIZE_P + N_ARRAY_SIZE_X + N_ARRAY_SIZE_Y , y_.begin());
        return y_;
    }

    void set_joint_para(joint_para_x_y_state_type joint_) {

        if (joint_.size() != (N_ARRAY_SIZE_P+N_ARRAY_SIZE_X+N_ARRAY_SIZE_Y)){
            std::cerr << "ERROR: set_joint_para: vector is not the right size" << std::endl;
        }

        //assign_concs();
        this->p = get_p_array(joint_);
        this->x = get_x_array(joint_);
        state_y_type y_ = get_y_array(joint_);
        set_ct_from_y(y_);
    }

    void set_joint_para_and_state(joint_para_x_y_state_type joint_) override{
        this->set_joint_para(joint_);
    }



    //! calculate all concentrations from amounts and conserved totals (ct vector)
    void assign_concs(){
#define ASSIGNMENT
#include MODEL_DEF_HEADER_FILE
#undef ASSIGNMENT
    }

    void initialise(){

        /*
        p = para_type();
        x_c = state_type();
        y_c = state_y_type();
        x = state_type();
        y = state_y_type();
        ct = CT_type();
        */

#define INITIAL
#include MODEL_DEF_HEADER_FILE
#undef INITIAL

#define FIXED
#include MODEL_DEF_HEADER_FILE
#undef FIXED

        assign_concs();


    }

    /*
    MODEL_CLASS_NAME(){
        initialise();

    }

    MODEL_CLASS_NAME(para_type p_init){
        initialise();
        p = p_init;
    }

    MODEL_CLASS_NAME(para_type p_init, CT_type ct_init){
        initialise();
        p = p_init;
        ct = ct_init;
    }
    */

    MODEL_CLASS_NAME(joint_para_x_y_state_type joint_) : p(N_ARRAY_SIZE_P), x_c(N_ARRAY_SIZE_X), x(N_ARRAY_SIZE_X), y_c(N_ARRAY_SIZE_Y), y(N_ARRAY_SIZE_Y), ct(N_ARRAY_SIZE_CT) {
        initialise();
        set_joint_para(joint_);
    }

    MODEL_CLASS_NAME() : p(N_ARRAY_SIZE_P), x_c(N_ARRAY_SIZE_X), x(N_ARRAY_SIZE_X), y_c(N_ARRAY_SIZE_Y), y(N_ARRAY_SIZE_Y), ct(N_ARRAY_SIZE_CT) {
        initialise();
        set_joint_para(new_joint_para_x_y_state_type());
    }

    virtual vector_type get_state() const override{
        return x;
    }


    state_y_type get_y_vector(const state_type &x_){
        x = x_;
        assign_concs();
        return y;
    }


    virtual vector_type get_dependent_vars(vector_type x_) override{
        return get_y_vector(x_);
    }


    state_y_type get_y_vector(const state_type &x_, const CT_type &ct_){
        ct = ct_;
        return get_y_vector(x_);
    }


    void set_ct_from_y(const state_y_type &y_){

        // clear ct[]
        CT_type::iterator it;
        for (it = ct.begin(); it != ct.end(); it++){
            *it = 0;
        }

        // calculate y when ct is zero
#define ASSIGNMENT
#include MODEL_DEF_HEADER_FILE
#undef ASSIGNMENT

        //calculate ct for y_
        for (unsigned int ii = 0 ; ii < ct.size(); ii++){
            ct[ii] = y_[ii] - y[ii];
        }

        // calculate y with new ct values
#define ASSIGNMENT
#include MODEL_DEF_HEADER_FILE
#undef ASSIGNMENT

        // verify values
        const double tol = 0.000001;
        for (unsigned int ii = 0 ; ii < y.size(); ii++){
            if (std::fabs(y_[ii] - y[ii]) > tol){
                std::cerr << "# !!! WARNING !!!: set_ct_from_y: y_[] and y[] vectors do not exactly match: diff["  << ii <<  "]=" << (y_[ii] - y[ii]) << std::endl;
            }
        }

    }

    static unsigned int get_p_index(const std::string name)  {
        if (name.compare("") == 0){ return NAME_NOT_FOUND; }
        std::vector<std::string>::const_iterator it = std::find(p_name_vec.begin(), p_name_vec.end(), name);
        if (it != p_name_vec.end()){
            return (it - p_name_vec.begin());
        }
        else {
            std::cerr << "ERROR: can't find parameter name:\t" << name << "\t in names[] vector" << std::endl;
            return NAME_NOT_FOUND;
        }
    }

    static unsigned int get_ct_index(const std::string name) {
        if (name.compare("") == 0){ return NAME_NOT_FOUND; }
        std::vector<std::string>::const_iterator it = std::find(ct_name_vec.begin(), ct_name_vec.end(), name);
        if (it != ct_name_vec.end()){
            return (it - ct_name_vec.begin());
        }
        else {
            std::cerr << "ERROR: can't find parameter name:\t" << name << "\t in ct_name_vec[] vector" << std::endl;
            return NAME_NOT_FOUND;
        }
    }

    static unsigned int get_x_index(const std::string name) {
        if (name.compare("") == 0){ return NAME_NOT_FOUND; }
        std::vector<std::string>::const_iterator it = std::find(x_name_vec.begin(), x_name_vec.end(), name);
        if (it != x_name_vec.end()){
            return (it - x_name_vec.begin());
        }
        else {
            std::cerr << "ERROR: can't find parameter name:\t" << name << "\t in x_name_vec[] vector" << std::endl;
            return NAME_NOT_FOUND;
        }
    }

    static unsigned int get_y_index(const std::string name) {
        if (name.compare("") == 0){ return NAME_NOT_FOUND; }
        std::vector<std::string>::const_iterator it = std::find(y_name_vec.begin(), y_name_vec.end(), name);
        if (it != y_name_vec.end()){
            return (it - y_name_vec.begin());
        }
        else {
            std::cerr << "ERROR: can't find parameter name:\t" << name << "\t in y_name_vec[] vector" << std::endl;
            return NAME_NOT_FOUND;
        }
    }

    static unsigned int get_joint_para_x_y_state_index(const std::string name) {
        if (name.compare("") == 0){ return NAME_NOT_FOUND; }

        std::vector<std::string>::const_iterator it = std::find(p_name_vec.begin(), p_name_vec.end(), name);
        if (it != p_name_vec.end()){
            return (it - p_name_vec.begin());
        }
        it = std::find(x_name_vec.begin(), x_name_vec.end(), name);
        if (it != x_name_vec.end()){
            return N_ARRAY_SIZE_P + (it - x_name_vec.begin());
        }
        it = std::find(y_name_vec.begin(), y_name_vec.end(), name);
        if (it != y_name_vec.end()){
            return N_ARRAY_SIZE_P + N_ARRAY_SIZE_X + (it - y_name_vec.begin());
        }

        std::cerr << "ERROR: can't find parameter name:\t" << name << "\t in p_names[], x_names[] or y_names[] vectors" << std::endl;
        return NAME_NOT_FOUND;
    }


    unsigned int get_joint_para_and_state_index(const std::string name) const override{
        return get_joint_para_x_y_state_index(name);
    }



    joint_para_x_y_state_type get_joint_para_and_state() override{
        assign_concs();
        return get_joint_para_x_y_state(this->p, this->x, this->y);
    }

    static unsigned int get_joint_x_y_state_index(const std::string name) {
        std::vector<std::string>::const_iterator it = std::find(x_name_vec.begin(), x_name_vec.end(), name);
        if (it != x_name_vec.end()){
            return (it - x_name_vec.begin());
        }
        it = std::find(y_name_vec.begin(), y_name_vec.end(), name);
        if (it != y_name_vec.end()){
            return N_ARRAY_SIZE_X + (it - y_name_vec.begin());
        }

        std::cerr << "ERROR: can't find parameter name:\t" << name << "\t in x_names[] or y_names[] vectors" << std::endl;
        return NAME_NOT_FOUND;
    }

    joint_x_y_state_type get_joint_x_y_state(){
        assign_concs();
        return get_joint_x_y_state(this->x, this->y);
    }

    virtual void initial_assignments() override {
    	assign_concs();

    	/*
    	for (auto it = y.begin() ; it != y.end(); it++){
    		std::cerr << *it << "\t";
    	}
    	std::cerr << std::endl;
    	*/

    	//INITIAL_ASSIGN
#define INITIAL_ASSIGN
#include MODEL_DEF_HEADER_FILE
#undef INITIAL_ASSIGN


    	/*
    	for (auto it = y.begin() ; it != y.end(); it++){
    		std::cerr << *it << "\t";
    	}
    	std::cerr << std::endl;
    	 */

    	// bug fix
    	state_y_type y_ = new_state_y_type();
    	y_ = y;
    	set_ct_from_y(y_);


    	/*
    	for (auto it = y.begin() ; it != y.end(); it++){
    		std::cerr << *it << "\t";
    	}
    	std::cerr << std::endl;
		*/


    	assign_concs();

    	/*
    	for (auto it = y.begin() ; it != y.end(); it++){
    		std::cerr << *it << "\t";
    	}
    	std::cerr << std::endl;
		*/


    }


    void operator()( const state_type &x_ , state_type &dx , double t ) override
    {
    x = x_;
    assign_concs();



    /*
    std::cout << "# t:\t" << t << std::endl;
    std::cout << "# x:";
    for (state_type::const_iterator it = x.begin(); it != x.end(); it++){
        std::cout << "\t" << *it;
    }
    std::cout << std::endl;

    std::cout << "# y:";
    for (state_y_type::const_iterator it = y.begin(); it != y.end(); it++){
        std::cout << "\t" << *it;
    }
    std::cout << std::endl;
    */


#define ODEs
#include MODEL_DEF_HEADER_FILE
#undef ODEs


    /*
    std::cout << "# dx:";
    for (state_type::const_iterator it = dx.begin(); it != dx.end(); it++){
        std::cout << "\t" << *it;
    }
    std::cout << std::endl;
    */


    }

    /*
    // this will not work as need to define functions as a pair<>() for odeint integration to work
    void first(const state_type &x_ , state_type &dx , double t ){
        return (*this)(x, dx, t);
    }

    void second(const state_type &x_ , state_matrix_type &jacobi, double t){
        const double h = 0.0001;
        jacobi = new_state_matrix_type();
        state_matrix_type low = new_state_matrix_type();
        state_matrix_type high = new_state_matrix_type();
    }
    */



};


#define NAME_ARRAYS
#include MODEL_DEF_HEADER_FILE
#undef NAME_ARRAYS
const std::vector<std::string> MODEL_CLASS_NAME::p_name_vec = std::vector<std::string>(p_names, carr_end(p_names));
const std::vector<std::string> MODEL_CLASS_NAME::x_name_vec = std::vector<std::string>(x_names, carr_end(x_names));
const std::vector<std::string> MODEL_CLASS_NAME::y_name_vec = std::vector<std::string>(y_names, carr_end(y_names));
const std::vector<std::string> MODEL_CLASS_NAME::xc_name_vec = std::vector<std::string>(xc_names, carr_end(xc_names));
const std::vector<std::string> MODEL_CLASS_NAME::pc_name_vec = std::vector<std::string>(pc_names, carr_end(pc_names));
const std::vector<std::string> MODEL_CLASS_NAME::yc_name_vec = std::vector<std::string>(yc_names, carr_end(yc_names));
const std::vector<std::string> MODEL_CLASS_NAME::dx_name_vec = std::vector<std::string>(dx_names, carr_end(dx_names));
const std::vector<std::string> MODEL_CLASS_NAME::ct_name_vec = std::vector<std::string>(ct_names, carr_end(ct_names));



}


#undef N_METABS
#undef N_ODE_METABS
#undef N_INDEP_METABS
#undef N_COMPARTMENTS
#undef N_GLOBAL_PARAMS
#undef N_KIN_PARAMS
#undef N_REACTIONS

#undef N_ARRAY_SIZE_P  	// number of parameters
#undef N_ARRAY_SIZE_X  	// number of initials
#undef N_ARRAY_SIZE_Y  	// number of assigned elements
#undef N_ARRAY_SIZE_XC 	// number of x concentration
#undef N_ARRAY_SIZE_PC 	// number of p concentration
#undef N_ARRAY_SIZE_YC 	// number of y concentration
#undef N_ARRAY_SIZE_DX 	// number of ODEs
#undef N_ARRAY_SIZE_CT 	// number of conserved totals

#undef T


//#endif /* BMEG_MODEL_H_ */

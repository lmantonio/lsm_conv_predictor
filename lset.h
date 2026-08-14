#ifndef LSET_H
#define LSET_H
#include "mesh.h"
#include "fem.h"
#include <string>

template <typename T> int sign(T val) {
    return (T(0) < val) - (val < T(0));
}

struct optInput{

	int hole_num_x;
	int hole_num_y;
	double hole_radius;
	int reset_lsf_freq;
	int reset_lsf_numStep;
	int solve_lsf_numStep_initial;
	int topologyOptimization_numIter_max;
	double lagrangeMultiplier_V;
	double penalizationParameter_V;
	double lagrangeMultiplier_P;
	double allowance_adv;
	double solve_lsf_top_freq;
	double solve_lsf_top_fracToRemove;
	double allowance_top;
	double volume_fraction_target;

};

//<new>
struct optMode{
	int optMode;
	std::string inputDensityFile;
};
//<new>

class lset{
	private:		
		int hole_numx;
    	int hole_numy;
		double hole_r;
		int num_step_lsf_reset;
		int reset_lsf_freq;
		int num_step_lsf_solve;
		int num_iter_max;
		double lagrangeMultiplier_V;
		double penalizationParameter_V;
		double lagrangeMultiplier_P;
		double allowance_adv;
		double allowance_top;
		int solve_lsf_top_freq;
		double frac_to_remove;
		double volFrac_target;
		int numNode;
		const nodeStruct *node;
		int numElement;
		const int *element;
		int nelx;
		int nely;
		double elementLengthX;
    	double elementLengthY;
		// <new>
		double *elementCompliance;
		double *elementCompliance_init;
		// <new>
		double *nodeCompliance;
		double *compliance;
		double *density;
		double lengthX;
    	double lengthY;		
		double ersatzMaterial;
		double CFL;
		double timeStep;		
		double volume;
		double volumeMax;
		double volumeTarget;
		double volumeFraction;
		double perimeter;
		double factorToAvoidZero;
		double factorToRegularizeNormalVelocity;
		// double objectiveFunction;
		// <new>
		double *objectiveFunction;
		// <new>
		double *lsf;
		double *lsf_mirror;		
		double *lsf_tmp;
		double *lsf_top;
		double *elementLsf;
		double *elementTriangleLsf;		
		double *normalVelocity;		
		double *normalVelocity_reg;
		double *normalVelocity_tmp;
		double *RHS_normalVelocity;
		double *lsf_mid_aux;
    	double *lsf_low_aux;
    	double *lsf_high_aux;				
		double *curvature;
		double *normalVector_x;
		double *normalVector_y;
		bool flag_topologicalGradient;		
		// <new>
		double *elementTopGrad;
		double *elementTopGrad_init;
		// <new>
		double *nodeTopGrad;
		bool seedHole;
		
		//variables for regularization matrix sparse assembly
		double *coo_val;
		int *coo_row;
		int *coo_col;
		int *coo_idx;
		double *csr_val;
		int *csr_col;
		int *csr_row;
		int csr_colsize;
		int csr_rowsize;
		double *cond_J;
		double *r;
		double *w;
		double *z;
		double *p;
		
		std::string dirPath;
		int simulationID;
		std::string meshFilePath;

		bool solverNotConverged;
		int sim_converged;

		//<new>
		int num_iteration;
		int lset_opt_mode_input;
		std::string initialTopology_file;
		bool flag_initialTopology;
		//<new>
		
		void set_lsf();
		void reset_lsf(int num_step_in_the_loop);
		void solve_lsf(int numStep);
		void solve_lsf_top();
		void set_regularizationMatrix_CSR();
		void solve_regularizationSystem_PCG();
		void set_normalVelocity();
		void regularize_normalVelocity();
		void set_normalVelocity_top();		
		void removeLsfZeros();
		double linearInterp();
		double set_dxm(int row, int col);
		double set_dxm_neumannGrad(int row, int col);
		double set_dxmxm(int row, int col);
		double set_dxmxm_neumannGrad(int row, int col);
		double set_dxp(int row, int col);
		double set_dxp_neumannGrad(int row, int col);
		double set_dxpxp(int row, int col);
		double set_dxpxp_neumannGrad(int row, int col);
		double set_dxmxp(int row, int col);
		double set_dxmxp_neumannGrad(int row, int col);
		double set_dym(int row, int col);
		double set_dym_neumannGrad(int row, int col);
		double set_dymym(int row, int col);
		double set_dymym_neumannGrad(int row, int col);
		double set_dyp(int row, int col);
		double set_dyp_neumannGrad(int row, int col);
		double set_dypyp(int row, int col);
		double set_dypyp_neumannGrad(int row, int col);
		double set_dymyp(int row, int col);
		double set_dymyp_neumannGrad(int row, int col);	
		double minmod(double arg1, double arg2);
		double flux(double u1,double u2, double v1, double v2);
		void set_timeStep();
		void set_curvature();
		void set_density();
		
		void set_volume();
		void set_volumeMax();
		void set_volumeFraction();
		void set_perimeter();
		void solve_lsf();
		double set_objectiveFunction();
		void set_compliance();
		void set_seedHole(int iterId);
		void set_lagrangeMultiplier_V();
		// bool check_topologyConvergence(int numIter, int lset_solutionAttempt);
		// <new>
		bool check_topologyConvergence(int numIter, int lset_solutionAttempt, int objectiveFunction_step);
		// <new>

		void log_volumeFraction();
		void log_objectiveFunction(double objectiveFunction);
		void log_iterationNumber(int iterNumber);
		void log_compliance();
		void log_density_hist();
		void log_lsf_hist();

		void log_data_finalresult(int iterNumber, bool convergenceOK, double simuTime);
		void log_simulationData(int iterNumber, bool convergenceOK, double simuTime);
		void log_density(std::string filename);
		void log_lsf(std::string filename);
		void remove_files(std::string filename);
		void log2files(std::string filename, int iter, double fob, double iter_compliance, double volfrac);
		void log2files_red(std::string filename, int iter, double fob, double iter_compliance, double volfrac);


	public:
		//lset(mesh &mesh2d, fem &femmodel, optInput lsetInput);
		//<new>
		lset(mesh &mesh2d, fem &femmodel, optInput lsetInput, optMode);
		//<new>
		~lset();				

		void optimize_topology(fem &femmodel);
		//<new>
		void get_initial_compliance_and_topgrad(fem &femmodel);
		//<new>
		inline int get_sim_converged() {return sim_converged;};
		inline double* get_density(){return density;};
		
		inline double* get_elementCompliance_init(){return elementCompliance_init;};	
		inline double* get_elementTopGrad_init(){return elementTopGrad_init;};	


		//<new>
		double get_compliance() {return compliance[0];};
		int get_num_iteration() {return num_iteration;};
		//<new>	
		
};

#endif
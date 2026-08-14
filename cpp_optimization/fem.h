#ifndef FEM_H
#define FEM_H
#include "mesh.h"
#include <string>

struct materialProp{
	double youngModulus;
	double poissonCoef;
};

class fem{
	private:
		int numNode;
		int numElement;
		int numBoundary;
		int numLoad;
		const nodeStruct *node;
		const int *element;
		const boundaryStruct *boundary;
		const loadStruct *load;
		double elementLengthX;
    	double elementLengthY;
		int nelx;
		int nely;		
		materialProp *material;
		double lameLambda;
		double lameMu;
		double *D;
		double *N;
		double *NdR;
		double *NdS;
		double *nodexCoord;
		double *nodeyCoord;
		double *J;
		double *invJ;
		double *B;
		double *tsme;
		double *ktr;
		double *B1;
		int *id;
		int eqNum;
		double *qe;
		double *u;
		int *loc;
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
		double* r;
		double* w;
		double* z;
		double* p;		
		double *ue;
		double *compliance;
		double *elementCompliance;
		double *nodeCompliance;		
		double *density;
		double *elementStrain;
		double *elementStress;
		double *elementTopGrad;
		double *nodeTopGrad;
		bool solverNotConverged;
		
		void matDef(int materialId);
		void set_materialProp(int materialId);
		void set_id_and_eqNum();
		void set_tsme();
		void set_density();		
		void assembly_CSR();
		void solver_pcg();

		void assembly_CSR_par();
		void solver_pcg_par();

		void set_qe(int loadId);
		void set_u();
		void set_ue(int elementId);
		void set_elementCompliance();
		void set_nodeCompliance();		
		void set_elementStrain(int elementId);
		void set_elementStress(int elementId);
		void set_nodeTopGrad();
		
		void log_node();
		void log_element();
		void log_boundary();
		void log_load();
		void log_tsme();
		void log_u();
		void log_u_nodeFormat();
		void log_elementCompliance();
		void log_nodeCompliance();
		void log_density();	
		
	public:
		fem(mesh &mesh2d);
		~fem();		
		
		void set_femAnalysis();
				
		inline double* get_compliance() {return compliance;};
		// <new>
		inline double* get_elementCompliance(){return elementCompliance;};
		// <new>
		inline double* get_nodeCompliance(){return nodeCompliance;};
		inline double* get_density(){return density;};
		// <new>
		inline double* get_elementTopGrad(){return elementTopGrad;};
		// <new>
		inline double* get_nodeTopGrad(){return nodeTopGrad;};
		inline bool get_solverNotConverged(){return solverNotConverged;}
					
};
#endif
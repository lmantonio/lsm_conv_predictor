#include "fem.h"
#include "mesh.h"
#include <cmath>
#include <iostream>
#include <algorithm>
#include <string>

#include <fstream>
#include <iomanip>
#include <sstream>

// ---windows only #include <omp.h>

#include <chrono>
#include <ctime>

//constructor for fem model
fem::fem(mesh &mesh2d):
	numNode(0),
	numElement(0),
	numBoundary(0),
	numLoad(0),
	node(nullptr),
	element(nullptr),
	boundary(nullptr),
	load(nullptr),
	elementLengthX(0),
    elementLengthY(0),
	nelx(0),
	nely(0),
	material(nullptr),
	lameLambda(0),
	lameMu(0),
	D(nullptr),
	N(nullptr),
	NdR(nullptr),
	NdS(nullptr),
	nodexCoord(nullptr),
	nodeyCoord(nullptr),
	J(nullptr),
	B(nullptr),
	tsme(nullptr),
	B1(nullptr),
	ktr(nullptr),
	id(nullptr),
	eqNum(0),
	qe(nullptr),
	u(nullptr),
	loc(nullptr),
	coo_val(nullptr),
	coo_row(nullptr),
	coo_col(nullptr),
	coo_idx(nullptr),
	csr_val(nullptr),
	csr_col(nullptr),
	csr_row(nullptr),
	cond_J(nullptr),
	r(nullptr),
	w(nullptr),
	z(nullptr),
	p(nullptr),	
	ue(nullptr),
	compliance(nullptr),
	elementCompliance(nullptr),
	nodeCompliance(nullptr),	
	density(nullptr),
	elementStrain(nullptr),
	elementStress(nullptr),
	elementTopGrad(nullptr),
	nodeTopGrad(nullptr),
	solverNotConverged(false){

	numNode = mesh2d.get_numNode();
	numElement = mesh2d.get_numElement();
	numBoundary = mesh2d.get_numBoundary();
	numLoad = mesh2d.get_numLoad();
	node = mesh2d.get_Nodes();
	element = mesh2d.get_Elements();
	boundary = mesh2d.get_Boundaries();
	load = mesh2d.get_Loads();
	elementLengthX = mesh2d.get_elementLengthX();
    elementLengthY = mesh2d.get_elementLengthY();
	nelx = mesh2d.get_nelx();
	nely = mesh2d.get_nely();

	material = new materialProp[1]();
	D = new double[9]();
	N = new double[4];
	NdR = new double[4];
	NdS = new double[4];
	nodexCoord = new double[4];
	nodeyCoord = new double[4];
	J = new double[4]();
	invJ = new double[4]();
	B = new double[24]();
	tsme = new double[64]();
	B1 = new double[8]();
	ktr = new double[64]();
	id = new int[numNode*2]();
	
	set_id_and_eqNum();

	u = new double[eqNum]();
	qe = new double[eqNum]();
	loc = new int[8]();
	coo_val = new double[64*numElement]();
	coo_row = new int[64*numElement]();
	coo_col = new int[64*numElement]();
	coo_idx = new int[64*numElement]();
	csr_val = new double[64*numElement]();
	csr_col = new int[64*numElement]();
	csr_row = new int[eqNum]();
	cond_J = new double[eqNum]();
	r = new double[eqNum]();
	w = new double[eqNum]();
	z = new double[eqNum]();
	p = new double[eqNum]();
	ue = new double[8]();
	compliance = new double[1]();
	elementCompliance = new double[numElement]();
	nodeCompliance = new double[numNode]();
	density = new double[numElement]();	
	elementStrain = new double[12]();
	elementStress = new double[12]();
	elementTopGrad = new double[numElement]();
	nodeTopGrad = new double[numNode]();
	
	set_materialProp(1);
	set_tsme();
	set_density();
	
}

//destructor for fem class
fem::~fem(){
	//delete [] node;
	//delete [] element;
	//delete [] boundary;
	//delete [] load;
	delete [] material;
	delete [] D;
	delete [] N;
	delete [] NdR;
	delete [] NdS;
	delete [] nodexCoord;
	delete [] nodeyCoord;
	delete [] J;
	delete [] B;
	delete [] tsme;
	delete [] B1;
	delete [] ktr;
	delete [] id;
	delete [] qe;
	delete [] u;
	delete [] loc;
	delete [] coo_val;
	delete [] coo_row;
	delete [] coo_col;
	delete [] coo_idx;
	delete [] csr_val;
	delete [] csr_col;
	delete [] csr_row;
	delete [] cond_J;
	delete [] r;
	delete [] w;
	delete [] z;
	delete [] p;
	delete [] ue;
	delete [] compliance;
	delete [] elementCompliance;
	delete [] nodeCompliance;
	delete [] density;	
	delete [] elementStress; 
	delete [] elementStrain;
	delete [] elementTopGrad;
	delete [] nodeTopGrad;
	
}

//=================================================================================================
// private member functions
//=================================================================================================

//material definition function (work as database)
void fem::matDef(int materialId){
	switch (materialId){
         case 1:
         	//steel 
          	material[0].youngModulus = 1;
          	material[0].poissonCoef = 0.3;
            break;
         case 2:
         	//second material ramdow chosen Kuow Bang test to validate
            material[0].youngModulus = 1e6;
          	material[0].poissonCoef = 0.3;
            break;
         default:
         	//steel 
          	material[0].youngModulus = 1;
          	material[0].poissonCoef = 0;            
      	}
}

void fem::set_materialProp(int materialId){
	matDef(materialId);
	lameLambda = material[0].youngModulus*material[0].poissonCoef/((1+material[0].poissonCoef)*(1-2*material[0].poissonCoef)); 
	lameMu = material[0].youngModulus/(2*(1+material[0].poissonCoef));

	D[0] = material[0].youngModulus/(1-material[0].poissonCoef*material[0].poissonCoef);
	D[1] = material[0].youngModulus*material[0].poissonCoef/(1-material[0].poissonCoef*material[0].poissonCoef);
	D[2] = 0;
	D[3] = D[1];
	D[4] = D[0];
	D[5] = 0;
	D[6] = 0;
	D[7] = 0;
	D[8] = 0.5*material[0].youngModulus/(1+material[0].poissonCoef);
	
}

void fem::set_id_and_eqNum(){
	for (int i=0 ; i<numBoundary; ++i){
		id[(boundary[i].nodeId-1)*2+(boundary[i].nodeDof-1)] = 1;
	}
	for (int i=0; i<numNode*2;++i){
		if (id[i] != 1){			
			id[i] = eqNum;
			++eqNum;								
		}else{
			id[i] = -1;
		}
	}
}

void fem::set_tsme(){
	//Gauss-Legendre quadrature for numerical integration
	double* gaussPt;
	double* gaussWt;
	double R;
	double S;
	double detJ;
	
	gaussPt = new double[8]();
	gaussWt = new double[8]();

	gaussPt[0] = -0.577350269189626;
	gaussPt[1] = -0.577350269189626;
	gaussPt[2] = +0.577350269189626;
	gaussPt[3] = -0.577350269189626;
	gaussPt[4] = +0.577350269189626;
	gaussPt[5] = +0.577350269189626;
	gaussPt[6] = -0.577350269189626;
	gaussPt[7] = +0.577350269189626;
	
	gaussWt[0] = 1;
	gaussWt[1] = 1;
	gaussWt[2] = 1;
	gaussWt[3] = 1;
	gaussWt[4] = 1;
	gaussWt[5] = 1;
	gaussWt[6] = 1;
	gaussWt[7] = 1;		
	
	for (int i=0; i<4; ++i){
		nodexCoord[i] = node[element[i+4]].xCoord;
		nodeyCoord[i] = node[element[i+4]].yCoord;
	}

	for (int i=0; i<4; ++i){
		R = gaussPt[i*2+0];
		S = gaussPt[i*2+1];

		N[0] = 0.25*(1-R)*(1-S);
		N[1] = 0.25*(1+R)*(1-S);
		N[2] = 0.25*(1+R)*(1+S);
		N[3] = 0.25*(1-R)*(1+S);
		NdR[0] = -0.25*(1-S);
		NdR[1] = 0.25*(1-S);
		NdR[2] = 0.25*(1+S);
		NdR[3] = -0.25*(1+S);
		NdS[0] = -0.25*(1-R);
		NdS[1] = -0.25*(1+R);
		NdS[2] = 0.25*(1+R);
		NdS[3] = 0.25*(1-R);

		J[0] = 0;
		J[1] = 0;
		J[2] = 0;
		J[3] = 0;

		for (int j=0; j<4; ++j){
			J[0] = J[0] + NdR[j]*nodexCoord[j];
			J[1] = J[1] + NdR[j]*nodeyCoord[j];
			J[2] = J[2] + NdS[j]*nodexCoord[j];
			J[3] = J[3] + NdS[j]*nodeyCoord[j];
		}

		detJ = J[0]*J[3]-J[2]*J[1];

		invJ[0] = (1/detJ)*J[3];
		invJ[1] = -(1/detJ)*J[1];
		invJ[2] = -(1/detJ)*J[2];
		invJ[3] = (1/detJ)*J[0];

		B[0] = invJ[0]*NdR[0]+invJ[1]*NdS[0];
		B[1] = 0;
		B[2] = invJ[0]*NdR[1]+invJ[1]*NdS[1];
		B[3] = 0;
		B[4] = invJ[0]*NdR[2]+invJ[1]*NdS[2];
		B[5] = 0;
		B[6] = invJ[0]*NdR[3]+invJ[1]*NdS[3];
		B[7] = 0;

		B[8] = 0;
		B[9] = invJ[2]*NdR[0]+invJ[3]*NdS[0];
		B[10] = 0;
		B[11] = invJ[2]*NdR[1]+invJ[3]*NdS[1];
		B[12] = 0;
		B[13] = invJ[2]*NdR[2]+invJ[3]*NdS[2];
		B[14] = 0;
		B[15] = invJ[2]*NdR[3]+invJ[3]*NdS[3];
		
		B[16] = B[9];
		B[17] = B[0];
		B[18] = B[11];
		B[19] = B[2];
		B[20] = B[13];
		B[21] = B[4];
		B[22] = B[15];
		B[23] = B[6];

		B1[0] = NdR[0]*invJ[0]+NdS[0]*invJ[1];
		B1[1] = NdR[0]*invJ[2]+NdS[0]*invJ[3];
		B1[2] = NdR[1]*invJ[0]+NdS[1]*invJ[1];
		B1[3] = NdR[1]*invJ[2]+NdS[1]*invJ[3];
		B1[4] = NdR[2]*invJ[0]+NdS[2]*invJ[1];
		B1[5] = NdR[2]*invJ[2]+NdS[2]*invJ[3];
		B1[6] = NdR[3]*invJ[0]+NdS[3]*invJ[1];
		B1[7] = NdR[3]*invJ[2]+NdS[3]*invJ[3];		
						
		//Element stiffness matrix construction | numerical integration
		double sum1(0);
		double sum2(0);

		for (int k=0; k<8; ++k){
			for (int l=0; l<8; ++l){
				for (int s=0; s<3; ++s){
					for (int r=0; r<3; ++r){
						sum1 = sum1 + D[s*3+r]*B[r*8+l];
					}
					sum2 = sum2 + sum1*B[s*8+k];
					sum1 = 0;
				}
				tsme[k*8+l] = tsme[k*8+l]+detJ*gaussWt[i*2+0]*gaussWt[i*2+1]*sum2;
				sum2 = 0;
				ktr[k*8+l] += detJ*gaussWt[i*2+0]*gaussWt[i*2+1]*(material[0].youngModulus/(1-material[0].poissonCoef))*B1[k]*B1[l];				
			}
		}		
	}

	delete [] gaussPt;
	delete [] gaussWt;
}

void fem::set_density(){
	for (int i = 0; i<numElement; ++i){
		density[i] = 1.0;
	}
}

void fem::assembly_CSR(){
	/*---------------------------------------------------------------
	Global FEM stiffness matrix assembly (TSM) in CSR (compress 
	sparse row) format.

	Algorithm strategy:
	01 - assembly TSM in COO format (coordinate sparse format);
		The COO format is composed by 3 arrays: 
		coo_val(double): array with the non-zeros values
		coo_col (int): array with col coordinates of non-zero values  
		coo_row (int): array with row coordinates of non-zero values

	02 - sort the 3 arrays driven by row order.
		This code uses std::sort() from library "algorithm", to sort
		a 4th index array. For this, the following class construction
		was done (from stack overflow):

		class sort_indices{
   			private:
     			int* mparr;
   			public:
     			sort_indices(int* parr) : mparr(parr) {}
     			bool operator()(int i, int j) const { return mparr[i]<mparr[j]; }
		};

		The 3 original arrays are sorted by the final index array.

	03 - sort the 3 arrays driven by column order, with the previus
		row ordering
		This code uses the same previus std::sort() function from 
		"algorithm" library to sort between columns segments (res-
		pecting the previus row ordering). 
		An auxiliary array is used in order to hold all the first new
		row entry. This will lead the segments first and last elments
		for the std::sort() function.
	
	04 - construct the CSR arrays by the following algorith:
		(i): check if there are entries for same row and columns.
			If this conditions is satisfied, make the summation of 
			all entries in the same situation.
		(ii): feed the corresponding column entry and feed the row array
		wherever a new row is observed with the value position.	
	---------------------------------------------------------------*/
	//COO format assembly
	int NZ(0);

	for (int i=0; i<this->numElement; ++i){		
		
		loc[0] = this->id[(this->element[i*8+4]-1)*2+0];
		loc[1] = this->id[(this->element[i*8+4]-1)*2+1];
		loc[2] = this->id[(this->element[i*8+5]-1)*2+0];
		loc[3] = this->id[(this->element[i*8+5]-1)*2+1];
		loc[4] = this->id[(this->element[i*8+6]-1)*2+0];
		loc[5] = this->id[(this->element[i*8+6]-1)*2+1];
		loc[6] = this->id[(this->element[i*8+7]-1)*2+0];
		loc[7] = this->id[(this->element[i*8+7]-1)*2+1];
			
		for (int m=0;m<8;++m){
			if (loc[m]>=0){
				for (int n=0;n<8;++n){
					if (loc[n]>=0){
						coo_val[NZ] = density[i]*tsme[m*8+n];						
						coo_row[NZ] = loc[n];
						coo_col[NZ] = loc[m]; 
						coo_idx[NZ] = NZ;
						++NZ;
					}
				}
			}
		}		
	}

	//sorting glogal stiffness matrix coo format by row, using std::sort and key sorting.
	//from stack overflow: for key sorting
	class sort_indices{
   		private:
    		int* mparr;
   		public:
     		sort_indices(int* parr) : mparr(parr) {}
     		bool operator()(int i, int j) const { return mparr[i]<mparr[j]; }
	};

	std::sort(coo_idx, coo_idx+NZ, sort_indices(coo_row));
	
	//#pragma omp parallel for num_threads(numThread)
	for (int i=0; i<NZ; ++i){
		csr_val[i] = coo_val[coo_idx[i]];
		csr_col[i] = coo_row[coo_idx[i]];
	}
	
	//#pragma omp parallel for num_threads(numThread)
	for (int i=0; i<NZ; ++i){
		coo_val[i] = csr_val[i];
		coo_row[i] = csr_col[i];
		csr_col[i] = coo_col[coo_idx[i]];
	}

	int counter = 0;
	csr_row[0] = 0;
	for (int i=0; i<NZ; ++i){
		coo_col[i] = csr_col[i];		
		if(coo_row[i] != coo_row[i+1]){
			++counter;
			csr_row[counter] = i+1;
		}
		coo_idx[i] = i;
	}
	++counter;
	csr_row[counter] = NZ;
	//-------------------------------------------------------------------->>>> implement parallel sort!!!
	for (int i=1; i<counter; ++i){
		std::sort(coo_idx+csr_row[i-1], coo_idx+csr_row[i], sort_indices(coo_col));		
	}
	//-------------------------------------------------------------------->>>> implement parallel sort!!!

	// ---windows only #pragma omp parallel for num_threads(numThread)
	for (int i=0; i<NZ; ++i){
		csr_val[i] = coo_val[coo_idx[i]];
		csr_col[i] = coo_col[coo_idx[i]];	
	}

	// ---windows only #pragma omp parallel for num_threads(numThread)
	for (int i=0; i<NZ; ++i){
		coo_val[i] = csr_val[i];
		coo_col[i] = csr_col[i];				
	}

	int NZV(0);
	counter = 0;

	csr_val[0] = coo_val[0];
	csr_col[0] = coo_col[0];
	csr_row[0] = 0;
	
	for (int i=1; i<NZ+1; ++i){
		if(coo_col[i] != coo_col[i-1] || coo_row[i] != coo_row[i-1]){			
			++NZV;
			csr_val[NZV] = coo_val[i];			
			csr_col[NZV] = coo_col[i];
			if(coo_row[i] != coo_row[i-1]){
				++counter;
				csr_row[counter] = NZV;				
			}			
		}else{
			csr_val[NZV] = csr_val[NZV] + coo_val[i];
		}
	}
	csr_colsize = NZV;
	csr_rowsize = counter;

	int k1; 
	int k2;

	counter = 0;	
	for (int i=0; i<csr_rowsize; ++i){
		k1 = csr_row[i];
		k2 = csr_row[i+1];
		for (int j=k1; j<k2;++j){
			 if(csr_col[j] == counter){
			 	cond_J[i] = 1.0/csr_val[j];
			 	++counter;
			 	break;
			 }
		}	
	}
	cond_J[csr_rowsize-1] = 1.0/csr_val[csr_col[csr_colsize]];
		
}

void fem::assembly_CSR_par(){
	/*---------------------------------------------------------------
	Global FEM stiffness matrix assembly (TSM) in CSR (compress 
	sparse row) format.

	Algorithm strategy:
	01 - assembly TSM in COO format (coordinate sparse format);
		The COO format is composed by 3 arrays: 
		coo_val(double): array with the non-zeros values
		coo_col (int): array with col coordinates of non-zero values  
		coo_row (int): array with row coordinates of non-zero values

	02 - sort the 3 arrays driven by row order.
		This code uses std::sort() from library "algorithm", to sort
		a 4th index array. For this, the following class construction
		was done (from stack overflow):

		class sort_indices{
   			private:
     			int* mparr;
   			public:
     			sort_indices(int* parr) : mparr(parr) {}
     			bool operator()(int i, int j) const { return mparr[i]<mparr[j]; }
		};

		The 3 original arrays are sorted by the final index array.

	03 - sort the 3 arrays driven by column order, with the previus
		row ordering
		This code uses the same previus std::sort() function from 
		"algorithm" library to sort between columns segments (res-
		pecting the previus row ordering). 
		An auxiliary array is used in order to hold all the first new
		row entry. This will lead the segments first and last elments
		for the std::sort() function.
	
	04 - construct the CSR arrays by the following algorith:
		(i): check if there are entries for same row and columns.
			If this conditions is satisfied, make the summation of 
			all entries in the same situation.
		(ii): feed the corresponding column entry and feed the row array
		wherever a new row is observed with the value position.	
	---------------------------------------------------------------*/
	//COO format assembly
	int NZ(0);

	for (int i=0; i<this->numElement; ++i){		
		
		loc[0] = this->id[(this->element[i*8+4]-1)*2+0];
		loc[1] = this->id[(this->element[i*8+4]-1)*2+1];
		loc[2] = this->id[(this->element[i*8+5]-1)*2+0];
		loc[3] = this->id[(this->element[i*8+5]-1)*2+1];
		loc[4] = this->id[(this->element[i*8+6]-1)*2+0];
		loc[5] = this->id[(this->element[i*8+6]-1)*2+1];
		loc[6] = this->id[(this->element[i*8+7]-1)*2+0];
		loc[7] = this->id[(this->element[i*8+7]-1)*2+1];
			
		for (int m=0;m<8;++m){
			if (loc[m]>=0){
				for (int n=0;n<8;++n){
					if (loc[n]>=0){
						coo_val[NZ] = density[i]*tsme[m*8+n];						
						coo_row[NZ] = loc[n];
						coo_col[NZ] = loc[m]; 
						coo_idx[NZ] = NZ;
						++NZ;
					}
				}
			}
		}		
	}

	//sorting glogal stiffness matrix coo format by row, using std::sort and key sorting.
	//from stack overflow: for key sorting
	class sort_indices{
   		private:
    		int* mparr;
   		public:
     		sort_indices(int* parr) : mparr(parr) {}
     		bool operator()(int i, int j) const { return mparr[i]<mparr[j]; }
	};

	std::sort(coo_idx, coo_idx+NZ, sort_indices(coo_row));
	
	// ---windows only #pragma omp parallel for //num_threads(4)
	for (int i=0; i<NZ; ++i){
		csr_val[i] = coo_val[coo_idx[i]];
		csr_col[i] = coo_row[coo_idx[i]];
	}
	
	// ---windows only #pragma omp parallel for //num_threads(4)
	for (int i=0; i<NZ; ++i){
		coo_val[i] = csr_val[i];
		coo_row[i] = csr_col[i];
		csr_col[i] = coo_col[coo_idx[i]];
	}

	int counter = 0;
	csr_row[0] = 0;
	for (int i=0; i<NZ; ++i){
		coo_col[i] = csr_col[i];		
		if(coo_row[i] != coo_row[i+1]){
			++counter;
			csr_row[counter] = i+1;
		}
		coo_idx[i] = i;
	}
	++counter;
	csr_row[counter] = NZ;
	//-------------------------------------------------------------------->>>> implement parallel sort!!!
	for (int i=1; i<counter; ++i){
		std::sort(coo_idx+csr_row[i-1], coo_idx+csr_row[i], sort_indices(coo_col));		
	}
	//-------------------------------------------------------------------->>>> implement parallel sort!!!

	// ---windows only #pragma omp parallel for //num_threads(4)
	for (int i=0; i<NZ; ++i){
		csr_val[i] = coo_val[coo_idx[i]];
		csr_col[i] = coo_col[coo_idx[i]];	
	}

	// ---windows only #pragma omp parallel for //num_threads(4)
	for (int i=0; i<NZ; ++i){
		coo_val[i] = csr_val[i];
		coo_col[i] = csr_col[i];				
	}

	int NZV(0);
	counter = 0;

	csr_val[0] = coo_val[0];
	csr_col[0] = coo_col[0];
	csr_row[0] = 0;
	
	for (int i=1; i<NZ+1; ++i){
		if(coo_col[i] != coo_col[i-1] || coo_row[i] != coo_row[i-1]){			
			++NZV;
			csr_val[NZV] = coo_val[i];			
			csr_col[NZV] = coo_col[i];
			if(coo_row[i] != coo_row[i-1]){
				++counter;
				csr_row[counter] = NZV;				
			}			
		}else{
			csr_val[NZV] = csr_val[NZV] + coo_val[i];
		}
	}
	csr_colsize = NZV;
	csr_rowsize = counter;

	int k1; 
	int k2;

	counter = 0;	
	for (int i=0; i<csr_rowsize; ++i){
		k1 = csr_row[i];
		k2 = csr_row[i+1];
		for (int j=k1; j<k2;++j){
			 if(csr_col[j] == counter){
			 	cond_J[i] = 1.0/csr_val[j];
			 	++counter;
			 	break;
			 }
		}	
	}
	cond_J[csr_rowsize-1] = 1.0/csr_val[csr_col[csr_colsize]];
		
}

void fem::solver_pcg(){
	/*---------------------------------------------------------------------------------------------
	Linear system solver using Left-preconditioned Conjugated Gradient
	Normal Residue (CGNR).	
	---------------------------------------------------------------------------------------------*/
	//matrix-vector multiplication
	int k1;
	int k2;
	int nummaxit(eqNum*2);//100//1000//eqNum
	double alp(0);
	double bet(0);
	double norm_w(0);
	double norm_b(0);
	double norm_r(0);
	double prod_zr1;
	double prod_zr2;	
	double prod1;
	double tol(0);
	double a;
	double b;
	double c;
	double delta;
	double ome;
	double aux1;
	double norm_u(0);
	double stopCriteria(1e-15);//1e-6//1e-10//1e-9//1e-12
	bool convergenceOK;

	convergenceOK = false;
	
	for (int i=0; i<csr_rowsize; ++i){
		
		k1 = csr_row[i];
		k2 = csr_row[i+1];
		prod1 = 0;

		for (int j=k1; j<k2;++j){
			prod1 += csr_val[j]*u[csr_col[j]];
		}

		w[i] = prod1;
		norm_b += qe[i]*qe[i];
		r[i] = qe[i] - prod1;
		z[i] = cond_J[i]*r[i];
		p[i] = z[i];
	}

	norm_b = sqrt(norm_b);
	prod_zr1 = 0;
	norm_r = 0;
	norm_w = 0;

	// ---windows only #pragma omp parallel for 
	for (int i=0; i<csr_rowsize; ++i){
		k1 = csr_row[i];
		k2 = csr_row[i+1];
		prod1=0;
		for (int j=k1; j<k2;++j){
			prod1 += csr_val[j]*p[csr_col[j]];
		}			
		w[i] = prod1;
		norm_w += prod1*p[i];
		norm_r += r[i]*r[i];
		prod_zr1 += z[i]*r[i];
	}		
	norm_r = sqrt(norm_r);
		
	alp = prod_zr1/norm_w;
	prod_zr2 = 0;
	norm_u = 0;
	
	// ---windows only #pragma omp parallel for 
	for (int i=0; i<csr_rowsize; ++i){
		u[i] += alp*p[i];
		r[i] -= alp*w[i]; 
		z[i] = cond_J[i]*r[i];
		prod_zr2 += z[i]*r[i];
		norm_u += u[i]*u[i];		
	}
	norm_u = sqrt(norm_u);
	bet = prod_zr2/prod_zr1;
	
	// ---windows only #pragma omp parallel for 
	for (int i=0; i<csr_rowsize; ++i){
		p[i] = z[i] + bet*p[i];
	}

	a = 1.0/alp;
	b = sqrt(bet/(alp*alp));
	c = 1.0;
	delta = a;

	for (int loop=1; loop < nummaxit; ++loop){
		prod_zr1 = 0;
		norm_r = 0;
		norm_w = 0;
		
		for (int i=0; i<csr_rowsize; ++i){
			k1 = csr_row[i];
			k2 = csr_row[i+1];
			prod1=0;
			for (int j=k1; j<k2;++j){
				prod1 += csr_val[j]*p[csr_col[j]];
			}			
			w[i] = prod1;
			norm_w += prod1*p[i];
			norm_r += r[i]*r[i];
			prod_zr1 += z[i]*r[i];
		}

		norm_r = sqrt(norm_r);
		aux1 = bet/alp;
		alp = prod_zr1/norm_w;
		a = 1.0/alp + aux1;
		ome = sqrt((delta-a)*(delta-a)+4.0*b*b*c*c);
		c = sqrt(0.5*(1-(delta-a)/ome));
		delta += ome*c*c;
		prod_zr2 = 0;
		norm_u = 0;
		
		// ---windows only #pragma omp parallel for
		for (int i=0; i<csr_rowsize; ++i){
			u[i] += alp*p[i];
			r[i] -= alp*w[i];
			z[i] = cond_J[i]*r[i];
			prod_zr2 += z[i]*r[i];
			norm_u += u[i]*u[i];
		}		
		
		norm_u = sqrt(norm_u);
		bet = prod_zr2/prod_zr1;
		b = sqrt(bet/(alp*alp));
		tol = stopCriteria*(delta*norm_u+norm_b);
		if (norm_r <= tol){
			convergenceOK = true;
			break;
		}
		//#pragma omp parallel for 
		for (int i=0; i<csr_rowsize; ++i){
			p[i] = z[i] + bet*p[i];
		}
	}

	if (!convergenceOK){
		std::cout << " - fem PCG Convergence ERROR";
		solverNotConverged = true;
	}
	
}

void fem::solver_pcg_par(){

	/*---------------------------------------------------------------------------------------------
	Linear system solver using Left-preconditioned Conjugated Gradient
	Normal Residue (CGNR).	
	---------------------------------------------------------------------------------------------*/
	//matrix-vector multiplication
	int k1;
	int k2;
	int nummaxit(eqNum*2);//100//1000//eqNum
	double alp(0);
	double bet(0);
	double norm_w(0);
	double norm_b(0);
	double norm_r(0);
	double prod_zr1;
	double prod_zr2;	
	double prod1;
	double tol(0);
	double a;
	double b;
	double c;
	double delta;
	double ome;
	double aux1;
	double norm_u(0);
	double stopCriteria(1e-15);//1e-6//1e-10//1e-9//1e-12 //--->Hugo: 1e-10
	bool convergenceOK;

	convergenceOK = false;
		 
	// ---windows only #pragma omp parallel for private(k1,k2,prod1) reduction(+:norm_b)
	for (int i=0; i<csr_rowsize; ++i){
		k1 = csr_row[i];
		k2 = csr_row[i+1];
		prod1 = 0;
		for (int j=k1; j<k2;++j){
			prod1 += csr_val[j]*u[csr_col[j]];
		}
		w[i] = prod1;
		norm_b += qe[i]*qe[i];
		r[i] = qe[i] - prod1;
		z[i] = cond_J[i]*r[i];
		p[i] = z[i];
	}
	
	norm_b = sqrt(norm_b);
	prod_zr1 = 0;
	norm_r = 0;
	norm_w = 0;
		
	// ---windows only #pragma omp parallel for private(k1,k2,prod1) reduction(+:norm_w) reduction(+:norm_r) reduction(+:prod_zr1) 
	for (int i=0; i<csr_rowsize; ++i){
		k1 = csr_row[i];
		k2 = csr_row[i+1];
		prod1=0;
		for (int j=k1; j<k2;++j){
			prod1 += csr_val[j]*p[csr_col[j]];
		}			
		w[i] = prod1;
		norm_w += prod1*p[i];
		norm_r += r[i]*r[i];
		prod_zr1 += z[i]*r[i];
	}		
		
	norm_r = sqrt(norm_r);		
	alp = prod_zr1/norm_w;
	prod_zr2 = 0;
	norm_u = 0;		
		
	// ---windows only #pragma omp parallel for private(alp) reduction(+:prod_zr2,norm_u)
	for (int i=0; i<csr_rowsize; ++i){
		u[i] += alp*p[i];
		r[i] -= alp*w[i]; 
		z[i] = cond_J[i]*r[i];
		prod_zr2 += z[i]*r[i];
		norm_u += u[i]*u[i];		
	}
		
	norm_u = sqrt(norm_u);
	bet = prod_zr2/prod_zr1;

	// ---windows only #pragma omp parallel for 
	for (int i=0; i<csr_rowsize; ++i){
		p[i] = z[i] + bet*p[i];
	}
	
	a = 1.0/alp;
	b = sqrt(bet/(alp*alp));
	c = 1.0;
	delta = a;

	for (int loop=1; loop < nummaxit; ++loop){
		prod_zr1 = 0;
		norm_r = 0;
		norm_w = 0;
				
		// ---windows only #pragma omp parallel for private(k1,k2,prod1) reduction(+:norm_w) reduction(+:norm_r) reduction(+:prod_zr1)
		for (int i=0; i<csr_rowsize; ++i){
			k1 = csr_row[i];
			k2 = csr_row[i+1];
			prod1=0;
			for (int j=k1; j<k2;++j){
				prod1 += csr_val[j]*p[csr_col[j]];
			}			
			w[i] = prod1;
			norm_w += prod1*p[i];
			norm_r += r[i]*r[i];
			prod_zr1 += z[i]*r[i];
		}

		norm_r = sqrt(norm_r);
		aux1 = bet/alp;
		alp = prod_zr1/norm_w;
		a = 1.0/alp + aux1;
		ome = sqrt((delta-a)*(delta-a)+4.0*b*b*c*c);
		c = sqrt(0.5*(1-(delta-a)/ome));
		delta += ome*c*c;
		prod_zr2 = 0;
		norm_u = 0;
	
		// ---windows only #pragma omp parallel for firstprivate(alp) reduction(+:prod_zr2,norm_u)
		for (int i=0; i<csr_rowsize; ++i){
			u[i] += alp*p[i];
			r[i] -= alp*w[i];
			z[i] = cond_J[i]*r[i];
			prod_zr2 += z[i]*r[i];
			norm_u += u[i]*u[i];
		}		
		
		norm_u = sqrt(norm_u);
		bet = prod_zr2/prod_zr1;
		b = sqrt(bet/(alp*alp));
		tol = stopCriteria*(delta*norm_u+norm_b);
		
		if (norm_r <= tol){
			convergenceOK = true;
			break;
		}

		// ---windows only #pragma omp parallel for
		for (int i=0; i<csr_rowsize; ++i){
			p[i] = z[i] + bet*p[i];
		}
	}

	if (!convergenceOK){
		std::cout << " - fem PCG Convergence ERROR";
		solverNotConverged = true;
	}

}

void fem::set_qe(int loadId){
	
	// ---windows only #pragma omp parallel for num_threads(4)
	for (int i=0; i<eqNum; ++i){
		qe[i] = 0;
	}

	if (id[(load[loadId].nodeId-1)*2+load[loadId].nodeDof-1] >= 0){
		qe[id[(load[loadId].nodeId-1)*2+load[loadId].nodeDof-1]] += load[loadId].loadModule;
	}	

}

void fem::set_u(){
	//assembly structure global stiffness matrix [K] - using compress sparse row
	//assembly_CSR();
	//solver linar system [K]{u} = {qe} to find {u}
		
	//solver_pcg();
	solver_pcg_par();

}

void fem::set_ue(int elementId){

	int dof1, dof2;
	int nodeToCheck;
	
	for (int j=0; j<4; ++j){

		nodeToCheck = element[(elementId-1)*8+(j+4)];
		
		dof1 = id[(nodeToCheck-1)*2];
		dof2 = id[(nodeToCheck-1)*2+1];

		if (dof1 != -1){
			ue[j*2] = u[dof1];
		}else{
			ue[j*2] = 0;
		}

		if (dof2 != -1){ 
			ue[j*2+1] = u[dof2];
		}else{
			ue[j*2+1] = 0;
		}

	}

}

void fem::set_elementCompliance(){
	double elemStrainEnergyDensity;
	double sum;
	
	for (int i = 0; i<numElement; ++i){
		set_ue(i+1);
		elemStrainEnergyDensity = 0;
		for (int j=0; j<8; ++j){
			sum = 0;
			for (int k=0; k<8; ++k){
				sum += ue[k]*tsme[8*j+k];	
			}
			elemStrainEnergyDensity += density[i]*sum*ue[j];
		}
		elementCompliance[i] += elemStrainEnergyDensity;

		compliance[0] += elemStrainEnergyDensity;
	}
}

void fem::set_nodeCompliance(){

    // zera vetor
    for (int i = 0; i < numNode; ++i){
        nodeCompliance[i] = 0.0;
    }

    // contador de quantos elementos contribuem por nó
    int* count = new int[numNode]();
    
    // loop nos elementos
    for (int e = 0; e < numElement; ++e){

        double comp = elementCompliance[e];

        // os 4 nós do elemento (últimas 4 posições)
        for (int j = 0; j < 4; ++j){

            int nodeId = element[e*8 + (4 + j)] - 1; // cuidado com indexação

            nodeCompliance[nodeId] += comp;
            count[nodeId] += 1;
        }
    }

    // média
    for (int i = 0; i < numNode; ++i){
        if (count[i] > 0){
            nodeCompliance[i] /= count[i];
        }
    }

    delete [] count;
}

/*
void fem::set_nodeCompliance(){

	double bottomLeft;
	double bottomRight;
	double topRight;
	double topLeft;
	
	// ---windows only #pragma omp parallel for private(bottomLeft,bottomRight,topRight,topLeft) num_threads(4)
	for (int row = 0; row < nely+1; ++row){
		for (int col = 0; col < nelx+1; ++col){
			
			bottomLeft = elementCompliance[(std::max(0,row-1))*nelx+std::max(0,col-1)];
			bottomRight = elementCompliance[(std::min(nely-1,row))*nelx+std::max(0,col-1)];
			topRight = elementCompliance[(std::min(nely-1,row))*nelx+std::min(nelx-1,col)];
			topLeft = elementCompliance[(std::max(0,row-1))*nelx+std::min(nelx-1,col)];

			nodeCompliance[row*(nelx+1)+col] = 0.25*(bottomLeft+bottomRight+topRight+topLeft);
		}
	}

}
*/	

void fem::set_elementStrain(int elementId){
	
	double R, S;
	double detJ;
	double sum;
	double *nodePts;

	nodePts = new double[8]();

	for (int i=0; i<4; ++i){
		nodexCoord[i] = node[element[(elementId-1)*8+(4+i)]].xCoord;
		nodeyCoord[i] = node[element[(elementId-1)*8+(4+i)]].yCoord;
	}
		
	nodePts[0] = 0.0;
	nodePts[1] = 0.0;
		
	for (int i=0; i<1; ++i){
		R = nodePts[i*2];
		S = nodePts[i*2+1];

		N[0] = 0.25*(1-R)*(1-S);
		N[1] = 0.25*(1+R)*(1-S);
		N[2] = 0.25*(1+R)*(1+S);
		N[3] = 0.25*(1-R)*(1+S);
		NdR[0] = -0.25*(1-S);
		NdR[1] = 0.25*(1-S);
		NdR[2] = 0.25*(1+S);
		NdR[3] = -0.25*(1+S);
		NdS[0] = -0.25*(1-R);
		NdS[1] = -0.25*(1+R);
		NdS[2] = 0.25*(1+R);
		NdS[3] = 0.25*(1-R);

		J[0] = 0;
		J[1] = 0;
		J[2] = 0;
		J[3] = 0;

		for (int j=0; j<4; ++j){
			J[0] = J[0] + NdR[j]*nodexCoord[j];
			J[1] = J[1] + NdR[j]*nodeyCoord[j];
			J[2] = J[2] + NdS[j]*nodexCoord[j];
			J[3] = J[3] + NdS[j]*nodeyCoord[j];
		}

		detJ = J[0]*J[3]-J[2]*J[1];

		invJ[0] = (1/detJ)*J[3];
		invJ[1] = -(1/detJ)*J[1];
		invJ[2] = -(1/detJ)*J[2];
		invJ[3] = (1/detJ)*J[0];

		B[0] = invJ[0]*NdR[0]+invJ[1]*NdS[0];
		B[1] = 0;
		B[2] = invJ[0]*NdR[1]+invJ[1]*NdS[1];
		B[3] = 0;
		B[4] = invJ[0]*NdR[2]+invJ[1]*NdS[2];
		B[5] = 0;
		B[6] = invJ[0]*NdR[3]+invJ[1]*NdS[3];
		B[7] = 0;
		B[8] = 0;
		B[9] = invJ[2]*NdR[0]+invJ[3]*NdS[0];
		B[10] = 0;
		B[11] = invJ[2]*NdR[1]+invJ[3]*NdS[1];
		B[12] = 0;
		B[13] = invJ[2]*NdR[2]+invJ[3]*NdS[2];
		B[14] = 0;
		B[15] = invJ[2]*NdR[3]+invJ[3]*NdS[3];
		B[16] = B[9];
		B[17] = B[0];
		B[18] = B[11];
		B[19] = B[2];
		B[20] = B[13];
		B[21] = B[4];
		B[22] = B[15];
		B[23] = B[6];

		for (int j=0; j<3; ++j){
			sum = 0;
			for (int k=0; k<8; ++k){
				sum += B[j*8+k]*ue[k];
			}
			elementStrain[i*3+j] = sum;						
		}

	}

	delete [] nodePts;
		
}

void fem::set_elementStress(int elementId){
	double sum;
	
	for (int i=0; i<1; ++i){	
		for (int j=0; j<3; ++j){
			sum = 0;
			for (int k=0; k<3; ++k){
				sum += D[j*3+k]*elementStrain[i*3+k];
			}
			elementStress[i*3+j] = density[elementId-1]*sum;
		}
	}
	
}

void fem::set_nodeTopGrad(){

    double aux1;
    double aux2;
    double aux3;

    aux1 = M_PI*(lameLambda+2*lameMu)/(2*lameMu)/(lameLambda+lameMu);

    // ===============================
    // 1. Calcula elementTopGrad
    // ===============================
    for (int i = 0; i < numElement; ++i){

        set_ue(i+1);
        set_elementStrain(i+1);
        set_elementStress(i+1);

        aux2 = 0;
        for (int j = 0; j < 3; ++j){
            aux2 += elementStrain[j]*elementStress[j];
        }

        aux2 = 4*lameMu*aux2;

        aux3 = (lameLambda-lameMu)*
               (elementStrain[0]+elementStrain[1])*
               (elementStress[0]+elementStress[1]);

        elementTopGrad[i] += aux1*(aux2 + aux3);
    }

    // ===============================
    // 2. Zera nodeTopGrad
    // ===============================
    for (int i = 0; i < numNode; ++i){
        nodeTopGrad[i] = 0.0;
    }

    // contador de contribuições
    int* count = new int[numNode]();

    // ===============================
    // 3. Espalha elemento → nós
    // ===============================
    for (int e = 0; e < numElement; ++e){

        double val = elementTopGrad[e];

        for (int j = 0; j < 4; ++j){

            int nodeId = element[e*8 + (4 + j)] - 1;

            nodeTopGrad[nodeId] += val;
            count[nodeId] += 1;
        }
    }

    // ===============================
    // 4. Média nos nós
    // ===============================
    for (int i = 0; i < numNode; ++i){
        if (count[i] > 0){
            nodeTopGrad[i] /= count[i];
        }
    }

    delete [] count;
}

/*
void fem::set_nodeTopGrad(){

	double aux1;
	double aux2;
	double aux3;
	
	aux1 = M_PI*(lameLambda+2*lameMu)/(2*lameMu)/(lameLambda+lameMu);

	for (int i=0; i<numElement;++i){

		set_ue(i+1);
		set_elementStrain(i+1);
		set_elementStress(i+1);

		aux2 = 0;
		for (int j=0; j<3; ++j){
			aux2 += elementStrain[j]*elementStress[j];
		}

		aux2 = 4*lameMu*aux2;

		aux3 = (lameLambda-lameMu)*(elementStrain[0]+elementStrain[1])*(elementStress[0]+elementStress[1]);

		elementTopGrad[i] += aux1*(aux2+aux3);		

	}

	double bottomLeft;
	double bottomRight;
	double topRight;
	double topLeft;

	for (int row = 0; row < nely+1; ++row){
		for (int col = 0; col < nelx+1; ++col){
			
			bottomLeft = elementTopGrad[(std::max(0,row-1))*nelx+std::max(0,col-1)];
			bottomRight = elementTopGrad[(std::min(nely-1,row))*nelx+std::max(0,col-1)];
			topRight = elementTopGrad[(std::min(nely-1,row))*nelx+std::min(nelx-1,col)];
			topLeft = elementTopGrad[(std::max(0,row-1))*nelx+std::min(nelx-1,col)];
	
			nodeTopGrad[row*(nelx+1)+col] = 0.25*(bottomLeft+bottomRight+topRight+topLeft);
		}
	}	
}
*/
//=================================================================================================
// public member functions
//=================================================================================================

void fem::set_femAnalysis(){
	
	compliance[0] = 0;
	for (int i=0; i<numElement; ++i){
		elementCompliance[i] = 0;		
		elementTopGrad[i] = 0;
	}	

	for (int i=0; i<numNode; ++i){
		nodeTopGrad[i] = 0;
	}

	//assembly_CSR();
	assembly_CSR_par();	
	
	for (int loadId=0; loadId<numLoad; ++loadId){
		if (solverNotConverged == false){
			set_qe(loadId);
			set_u();		
			set_nodeTopGrad();
			set_elementCompliance();
		}
	}
	set_nodeCompliance();			

}

void fem::log_node(){
	std::ofstream csv_file_nodes;
	csv_file_nodes.open("result-repo/nodes_fem.csv");

	csv_file_nodes << "nodeId" << "," << "x" << "," << "y" << "\n";
	for (int row = 0; row < numNode; ++row){
		csv_file_nodes << node[row].nodeId << ",";
		csv_file_nodes << node[row].xCoord << ",";
		csv_file_nodes << node[row].yCoord << "\n";
	}
	csv_file_nodes.close();
}

void fem::log_element(){
	std::ofstream csv_file_elements;
	csv_file_elements.open("result-repo/elements_fem.csv");

	for (int row = 0; row < numElement; ++row){
		for (int col = 0; col < 7; ++col){
			csv_file_elements << element[row*8+col] << ",";
		}
		csv_file_elements << element[row*8+7] << "\n";
	}
	csv_file_elements.close();
}

void fem::log_boundary(){
	std::ofstream csv_file_boundaries;
	csv_file_boundaries.open("result-repo/boundaries_fem.csv");

	csv_file_boundaries << "nodeId" << "," << "nodeDof" << "," << "boundaryModule" << "\n";
	for (int row = 0; row < numBoundary; ++row){
		csv_file_boundaries << boundary[row].nodeId << ",";
		csv_file_boundaries << boundary[row].nodeDof << ",";
		csv_file_boundaries << boundary[row].boundaryModule << "\n";
	}
	csv_file_boundaries.close();
}

void fem::log_load(){
	std::ofstream csv_file_loads;
	csv_file_loads.open("result-repo/loads_fem.csv");

	csv_file_loads << "nodeId" << "," << "nodeDof" << "," << "loadModule" << "\n";
	for (int row = 0; row < numLoad; ++row){
		csv_file_loads << load[row].nodeId << ",";
		csv_file_loads << load[row].nodeDof << ",";
		csv_file_loads << load[row].loadModule << "\n";
	}
	csv_file_loads.close();
}

void fem::log_tsme(){
	std::ofstream csv_file_tsme;
	csv_file_tsme.open("result-repo/tsme.csv");

	for (int row = 0; row < 8; ++row){
		for (int col = 0; col < 7; ++col){
			csv_file_tsme << tsme[row*8+col] << ",";
		}
		
		csv_file_tsme << tsme[row*8+7] << "\n";
	}

	csv_file_tsme.close();
}

void fem::log_u(){

	std::ofstream myufile;
	myufile.open ("result-repo/u.csv");

	for (int row = 0; row < eqNum; ++row){
		myufile << u[row] << std::endl;
	}

	myufile.close();

}

void fem::log_u_nodeFormat(){
	std::ofstream myufile;
	myufile.open ("result-repo/u_nodeFormat.csv");
	for (int i=0; i<numNode; ++i){

		if (id[i*2] >= 0){
			myufile << std::setprecision(10) << u[id[i*2]] << ",";
		}else{
			myufile << std::setprecision(10) << 0 << ",";
		}

		if (id[i*2+1] >= 0){
			myufile << std::setprecision(10) << u[id[i*2+1]] << std::endl;
		}else{
			myufile << std::setprecision(10) << 0 << std::endl;
		}
		
	}

	myufile.close();
}

void fem::log_elementCompliance(){
	std::ofstream myelementcompliancefile;
	myelementcompliancefile.open ("result-repo/elementCompliance.csv");
	for (int i=0; i<numElement; ++i){		
		myelementcompliancefile << std::setprecision(10) << elementCompliance[i] << std::endl;				
	}

	myelementcompliancefile.close();
}

void fem::log_nodeCompliance(){
	std::ofstream mynodecompliancefile;
	mynodecompliancefile.open ("result-repo/nodeCompliance.csv");
	for (int i=0; i<numNode; ++i){		
		mynodecompliancefile << std::setprecision(18) << nodeCompliance[i] << std::endl;				
	}

	mynodecompliancefile.close();
}

void fem::log_density(){
	std::ofstream mynodecompliancefile;
	mynodecompliancefile.open ("result-repo/density-fem.csv");
	for (int i=0; i<numElement; ++i){		
		mynodecompliancefile << std::setprecision(10) << density[i] << std::endl;				
	}

	mynodecompliancefile.close();
}


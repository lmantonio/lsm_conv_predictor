#include "lset.h"
#include "fem.h"
#include "mesh.h"
#include <cmath>
#include <iostream>
#include <algorithm>
#include <chrono> 
#include <ctime>
#include <string>

#include <fstream>
#include <iomanip>
#include <sstream>

#include <chrono>
#include <ctime>

//constructor for level set topology optimization (shape gradient)
//lset::lset(mesh &mesh2d, fem &femmodel, optInput lsetInput:
//<new>
lset::lset(mesh &mesh2d, fem &femmodel, optInput lsetInput, optMode lsetMode):
//<new>
    hole_numx(0),
    hole_numy(0),
    hole_r(0),
    num_step_lsf_reset(0),
    reset_lsf_freq(0),
    num_step_lsf_solve(0),    
    num_iter_max(0),
    lagrangeMultiplier_V(0),
    penalizationParameter_V(0),
    lagrangeMultiplier_P(0),
    allowance_adv(0),    
    allowance_top(0),
    solve_lsf_top_freq(0),
    frac_to_remove(0),
    volFrac_target(0),
    numNode(0),
    node(nullptr),
    numElement(0),
    element(nullptr),
    nelx(0),
    nely(0),
    elementLengthX(0),
    elementLengthY(0),
    // <new>
    elementCompliance(nullptr),
    elementCompliance_init(nullptr),
    // <new>    
    nodeCompliance(nullptr),
    compliance(nullptr),
    density(nullptr),
    // <new>
    elementTopGrad(nullptr),
    elementTopGrad_init(nullptr),
    // <new>
    nodeTopGrad(nullptr),   
    lengthX(0),
    lengthY(0),
    ersatzMaterial(1e-3),   
    CFL(1),
    timeStep(0),
    volume(0),
    volumeMax(0),
    volumeTarget(0),
    volumeFraction(0),
    perimeter(0),
    factorToAvoidZero(0),
    factorToRegularizeNormalVelocity(0),
    // objectiveFunction(0),
    // <new>
    objectiveFunction(nullptr),
    // <new>
    lsf(nullptr),
    lsf_mirror(nullptr),
    lsf_tmp(nullptr),
    lsf_top(nullptr),
    elementLsf(nullptr),
    elementTriangleLsf(nullptr),
    normalVelocity(nullptr),
    normalVelocity_reg(nullptr),
    normalVelocity_tmp(nullptr),
    RHS_normalVelocity(nullptr),
    lsf_mid_aux(nullptr),
    lsf_low_aux(nullptr),
    lsf_high_aux(nullptr),    
    curvature(nullptr),
    normalVector_x(nullptr),
    normalVector_y(nullptr),
    flag_topologicalGradient(true),
    seedHole(false),    
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
    dirPath("no path defined"),
    simulationID(0),
    meshFilePath("no mesh path def"),
    sim_converged(0),
    solverNotConverged(false),
    // <new>
    num_iteration(0),
    lset_opt_mode_input(0),
    initialTopology_file("None"),
    flag_initialTopology(false){
    // <new>
    //level set inputs
    hole_numx = lsetInput.hole_num_x;
    hole_numy = lsetInput.hole_num_y;
    hole_r = lsetInput.hole_radius;
    num_step_lsf_reset = lsetInput.reset_lsf_numStep;
    reset_lsf_freq = lsetInput.reset_lsf_freq;
    num_step_lsf_solve = lsetInput.solve_lsf_numStep_initial;   
    num_iter_max = lsetInput.topologyOptimization_numIter_max;
    lagrangeMultiplier_V = lsetInput.lagrangeMultiplier_V;
    penalizationParameter_V = lsetInput.penalizationParameter_V;
    lagrangeMultiplier_P = lsetInput.lagrangeMultiplier_P;
    allowance_adv = lsetInput.allowance_adv;
    allowance_top = lsetInput.allowance_top;    
    solve_lsf_top_freq = lsetInput.solve_lsf_top_freq;
    frac_to_remove = lsetInput.solve_lsf_top_fracToRemove;
    volFrac_target = lsetInput.volume_fraction_target;
    
    // <new>
    lset_opt_mode_input = lsetMode.optMode;
    initialTopology_file = lsetMode.inputDensityFile;
    // <new>

    numNode = mesh2d.get_numNode();
    node = mesh2d.get_Nodes();
    numElement = mesh2d.get_numElement();
    element = mesh2d.get_Elements();
    nelx = mesh2d.get_nelx();
    nely = mesh2d.get_nely();
    elementLengthX = mesh2d.get_elementLengthX();
    elementLengthY = mesh2d.get_elementLengthY();
    
    // <new>
	elementCompliance = femmodel.get_elementCompliance();
    elementCompliance_init = new double[numElement]();
    for (int i=0; i<(numElement);++i){
        elementCompliance_init[i] = 0;
    }
	// <new>

    nodeCompliance = femmodel.get_nodeCompliance();
    compliance = femmodel.get_compliance();
    density = femmodel.get_density();
    // <new>
	elementTopGrad = femmodel.get_elementTopGrad();
    elementTopGrad_init = new double[numElement]();
    for (int i=0; i<(numElement);++i){
        elementTopGrad_init[i] = 0;
    }
	// <new>
    nodeTopGrad = femmodel.get_nodeTopGrad();
    solverNotConverged = femmodel.get_solverNotConverged();    

    lengthX = static_cast <double>(nelx)*elementLengthX;
    lengthY = static_cast <double>(nely)*elementLengthY;
    
    set_volumeMax();
    set_volume();
    set_volumeFraction();

    factorToAvoidZero = std::min(elementLengthX,elementLengthY)/20;
    factorToRegularizeNormalVelocity = 4*pow(elementLengthX,2);

    // <new>
    objectiveFunction = new double[num_iter_max+20]();
    // <new>

    lsf = new double[numNode]();
    lsf_mirror = new double[numNode]();
    lsf_tmp = new double[numNode]();
    lsf_top = new double[numNode]();
    elementLsf = new double[5]();
    elementTriangleLsf = new double[12]();
    normalVelocity = new double[numNode]();
    normalVelocity_reg = new double[numNode]();
    normalVelocity_tmp = new double[numNode]();
    RHS_normalVelocity = new double[numNode]();
    lsf_mid_aux = new double [3]();
    lsf_low_aux = new double [3]();
    lsf_high_aux = new double [3]();
    curvature = new double[numNode]();
    normalVector_x = new double[numNode]();
    normalVector_y = new double[numNode]();
    coo_val = new double[numNode*6]();
	coo_row = new int[numNode*6]();
	coo_col = new int[numNode*6]();
	coo_idx = new int[numNode*6]();
	csr_val = new double[numNode*6]();
	csr_col = new int[numNode*6]();
	csr_row = new int[numNode*6](); 
    cond_J = new double[numNode];
    r = new double[numNode]();
	w = new double[numNode]();
	z = new double[numNode]();
	p = new double[numNode]();    

    //Initialize density
    // if (solve_lsf_top_freq == 0) {
    //     flag_topologicalGradient = false;
    // }

    //<new>
    if (lset_opt_mode_input == 0){
        flag_topologicalGradient = false;
        solve_lsf_top_freq = 0;
    }else{
        if (solve_lsf_top_freq == 0){
            flag_topologicalGradient = false;
        }else{
            flag_topologicalGradient = true;
        }
    }
    //<new>

    // função para recuperar lset nos nós a partir da densidade do elemento
    set_lsf();

    set_regularizationMatrix_CSR();
    
    meshFilePath = mesh2d.get_meshFilePath();
    
}

//destructor for level set topology optimization
lset::~lset(){
	  	
	// delete [] node;
    // delete [] element;	
    // delete [] nodeCompliance;
    // delete [] compliance;
    // delete [] density;
    // delete [] nodeTopGrad;
	delete [] lsf;
    delete [] lsf_mirror;
    delete [] lsf_tmp;
    delete [] lsf_top;
    delete [] elementLsf;
    delete [] elementTriangleLsf;
    delete [] normalVelocity;
    delete [] normalVelocity_reg;
    delete [] normalVelocity_tmp;
	delete [] RHS_normalVelocity;
    delete [] lsf_mid_aux;
    delete [] lsf_low_aux;
    delete [] lsf_high_aux;
    delete [] curvature;
    delete [] normalVector_x;
    delete [] normalVector_y;
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
    delete [] elementCompliance_init;
    delete [] elementTopGrad_init;
    
}

//=================================================================================================
// Private member functions
//=================================================================================================

void lset::set_lsf(){

    // if (flag_topologicalGradient == false){

    //     double auxY, auxX, temp_lsf;

    //     auxY = ((hole_numy+1)*M_PI)/lengthY;
    //     auxX = ((hole_numx+1)*M_PI)/lengthX;
       
    //     for (int i=0; i<numNode; ++i){        
    //         temp_lsf = -cos(auxY*node[i].yCoord)*cos(auxX*node[i].xCoord)+hole_r-1;        
    //         if (fabs(temp_lsf)<1e-10){
    //             temp_lsf = 0.0;            
    //         }        
    //         lsf[i] = temp_lsf;
    //         lsf[i] = 0.2*ceil(std::max(temp_lsf,0.0))-0.1;        
    //     }

    // }else{
    //     for (int i=0; i<numNode; ++i){
    //         lsf[i] = -1.0;
    //     }
    // }

    //<new>
    std::ifstream input_top_file;
    input_top_file.open(initialTopology_file, std::ios::in);
    
    if (input_top_file){
        flag_initialTopology = true;
        //03 - based on element density, estimate value 
        std::string str_tmp;
        double ALP;
        double BET;
        double density_average;
        
        for (int i=0; i<(numElement-1); ++i){
            std::getline(input_top_file, str_tmp,',');
            density[i] = stold(str_tmp);
        }
        std::getline(input_top_file, str_tmp);
        density[numElement-1] = stold(str_tmp);

        ALP = -2.0020020020020022; 
        BET = +1.0020020020020022;

        for (int i=0; i<numNode; ++i){
            density_average = 0;
            for (int j=0; j<(numElement); ++j){
                if (element[j*8+4]==(i+1) ||
                    element[j*8+5]==(i+1) ||
                    element[j*8+6]==(i+1) ||
                    element[j*8+7]==(i+1)){
                        density_average += density[j];
                    }
            }
            density_average = 0.25*density_average;

            lsf[i] = ALP*density_average + BET;
        }

    }else{
        flag_initialTopology = false;
        if (flag_topologicalGradient == false){

            double auxY, auxX, temp_lsf;

            auxY = ((hole_numy+1)*M_PI)/lengthY;
            auxX = ((hole_numx+1)*M_PI)/lengthX;
       
            for (int i=0; i<numNode; ++i){        
                temp_lsf = -cos(auxY*node[i].yCoord)*cos(auxX*node[i].xCoord)+hole_r-1;        
                if (fabs(temp_lsf)<1e-10){
                    temp_lsf = 0.0;            
                }        
                lsf[i] = temp_lsf;
                lsf[i] = 0.2*ceil(std::max(temp_lsf,0.0))-0.1;        
            }

        }else{
            for (int i=0; i<numNode; ++i){
                lsf[i] = -1.0;
            }
        }
    }
    input_top_file.close();
    //<new>
        
    reset_lsf(50); // hard coded as 50, to reduce the number of users inputs
    set_density();
    
}

void lset::reset_lsf(int num_step_in_the_loop){
    //function to time march level set function - reinitializing
    double dxm, dxp, dxmxm, dxpxp, dxmxp;
    double dym, dyp, dymym, dypyp, dymyp;    

    double flux_aux_A, flux_aux_B, flux_aux_C, flux_aux_D;
    double delp2, delm2;
    double nabla, lsf_sign;
    double ds,dt;

    dt = 0.5*(std::min(elementLengthX,elementLengthY));
    ds = 0.05*sqrt(pow(elementLengthX,2)+pow(elementLengthY,2));
        
    for (int n=0; n<num_step_in_the_loop; ++n){

       
        // ---windows only #pragma omp parallel for
        for (int i=0; i<numNode; ++i){
            lsf_mirror[i] = lsf[i];
        }        
        
        // ---windows only #pragma omp parallel for private(dxm,dxp,dym,dyp,dxmxm,dxpxp,dxmxp,dymym,dypyp,dymyp,flux_aux_A,flux_aux_B,flux_aux_C,flux_aux_D,delp2,delm2,nabla,lsf_sign) collapse(2)
        for (int row=0; row<(nely+1); ++row){            
            for (int col=0; col<(nelx+1); ++col){

                dxm = set_dxm(row,col);              
                dxp = set_dxp(row,col);
                dym = set_dym(row,col);
                dyp = set_dyp(row,col);
                dxmxm = set_dxmxm(row,col);
                dxpxp = set_dxpxp(row,col);
                dxmxp = set_dxmxp(row,col);
                dymym = set_dymym(row,col);
                dypyp = set_dypyp(row,col);
                dymyp = set_dymyp(row,col);
                
                flux_aux_A = dxm + 0.5*elementLengthX*minmod(dxmxm,dxmxp);
                flux_aux_B = dxp - 0.5*elementLengthX*minmod(dxpxp,dxmxp);
                flux_aux_C = dym + 0.5*elementLengthY*minmod(dymym,dymyp);
                flux_aux_D = dyp - 0.5*elementLengthY*minmod(dypyp,dymyp);

                delp2 = flux(flux_aux_A,flux_aux_B,flux_aux_C,flux_aux_D);
                delm2 = flux(flux_aux_B,flux_aux_A,flux_aux_D,flux_aux_C);                     

                nabla = (pow(dxm,2)+pow(dxp,2)+pow(dym,2)+pow(dyp,2))*ds;

                lsf_sign = lsf_mirror[row*(nelx+1)+col]/sqrt(pow(lsf_mirror[row*(nelx+1)+col],2)+nabla);
                                
                lsf[row*(nelx+1)+col] += -dt*( std::max(lsf_sign,0.0)*delp2 + std::min(lsf_sign,0.0)*delm2-lsf_sign );         
                
            }            
        }//loop de varredura da matriz

    }//loop de evolução da lsf

}

void lset::solve_lsf(int numStep){

    double dxm, dxp, dxmxm, dxpxp, dxmxp;
    double dym, dyp, dymym, dypyp, dymyp;    

    double flux_aux_A, flux_aux_B, flux_aux_C, flux_aux_D;
    double delp2, delm2;

    double Vm, Vp;
    double gradx_lsf, grady_lsf, mag;
    
    for (int i=0; i<numStep; ++i){
        
        if ((i+1)%reset_lsf_freq==0){
            reset_lsf(num_step_lsf_reset);
        }
        
        set_curvature();

        // ---windows only #pragma omp parallel for
        for(int i=0; i<numNode; ++i){
            lsf_mirror[i] = lsf[i];
        }

        // ---windows only #pragma omp parallel for private(dxm,dxp,dym,dyp,dxmxm,dxpxp,dxmxp,dymym,dypyp,dymyp,flux_aux_A,flux_aux_B,flux_aux_C,flux_aux_D,delp2,delm2,gradx_lsf,grady_lsf,mag,Vp,Vm) collapse(2)
        for (int row=0; row<(nely+1); ++row){
            for (int col=0; col<(nelx+1); ++col){
        
                dxm = set_dxm_neumannGrad(row,col);              
                dxp = set_dxp_neumannGrad(row,col);
                dym = set_dym_neumannGrad(row,col);
                dyp = set_dyp_neumannGrad(row,col);                

                dxmxm = set_dxmxm_neumannGrad(row,col);
                dxpxp = set_dxpxp_neumannGrad(row,col);
                dxmxp = set_dxmxp_neumannGrad(row,col);
                dymym = set_dymym_neumannGrad(row,col);
                dypyp = set_dypyp_neumannGrad(row,col);
                dymyp = set_dymyp_neumannGrad(row,col);
                
                flux_aux_A = dxm + 0.5*elementLengthX*minmod(dxmxm,dxmxp);
                flux_aux_B = dxp - 0.5*elementLengthX*minmod(dxpxp,dxmxp);
                flux_aux_C = dym + 0.5*elementLengthY*minmod(dymym,dymyp);
                flux_aux_D = dyp - 0.5*elementLengthY*minmod(dypyp,dymyp);

                delp2 = flux(flux_aux_A,flux_aux_B,flux_aux_C,flux_aux_D);
                delm2 = flux(flux_aux_B,flux_aux_A,flux_aux_D,flux_aux_C);                     

                gradx_lsf = 0.5*(set_dxp_neumannGrad(row,col)+set_dxm_neumannGrad(row,col));
                grady_lsf = 0.5*(set_dyp_neumannGrad(row,col)+set_dym_neumannGrad(row,col));

                mag = sqrt(pow(gradx_lsf,2)+pow(grady_lsf,2)+pow(factorToAvoidZero,2)); 

                Vp = std::max(normalVelocity[row*(nelx+1)+col],0.0);
                Vm = std::min(normalVelocity[row*(nelx+1)+col],0.0);
                
                lsf[row*(nelx+1)+col] += -timeStep*(delp2*Vp+delm2*Vm) + timeStep*lagrangeMultiplier_P*curvature[row*(nelx+1)+col]*mag;                

            }
        }//loop de varredura da matriz

    }
    
}

void lset::solve_lsf_top(){
    
    double min_velocity(1e30);
    double perv;
    double perv_max;
    double perv_min;
    double volume_in;
    double volume_target;
    double volume_test;
    double counter(0);

    set_volume();
    volume_in = volume;

    for (int i=0; i<numNode; ++i){
        if (lsf[i] < 0 && normalVelocity[i] < 0){
            min_velocity = std::min(min_velocity,normalVelocity[i]);
        }
        lsf_top[i] = lsf[i];        
    }

    if (min_velocity != 0){

        perv = 0.95;
        volume_target = (1-frac_to_remove)*volume_in;
        volume_test = 0;     
        counter = 0;

        for (int i=0; i<numNode; ++i){
            if ( (lsf_top[i] < 0) && (normalVelocity[i] < 0) && (fabs(normalVelocity[i]) > perv*fabs(min_velocity)) ) {
                lsf[i] = -lsf_top[i];                
            }       
        }

        reset_lsf(100);
        set_density();
        set_volume();
        volume_test = volume;

        if ( (fabs(volume_test-volume_target)/volume_in) > 0.005 ){
            if(volume_test > volume_target){

                perv_max = perv;
                perv_min = perv;

                while ( (volume_test > volume_target) && (perv_min > 0) ){
                    perv_max = perv_min;
                    perv_min = std::max(perv_min-0.01,0.0);
                    for (int i=0; i<numNode; ++i){
                        if ( (lsf_top[i]<0) && (normalVelocity[i]<0) && (fabs(normalVelocity[i]) > perv_min*fabs(min_velocity)) ){
                            lsf[i] = -lsf_top[i];
                        }else{
                            lsf[i] = lsf_top[i];
                        }
                    }
                    reset_lsf(100);
                    set_density();
                    set_volume();
                    volume_test = volume;
                }

            }else{

                perv_max = perv;
                perv_min = perv;
                
                while( (volume_test <= volume_target) && (perv_max < 1.0) ){
                    perv_min = perv_max;
                    perv_max = std::min(perv_max+0.01,1.0);
                    for(int i=0; i<numNode; ++i){
                        if ( (lsf_top[i]<0) && (normalVelocity[i]<0) && (fabs(normalVelocity[i]) > perv_max*fabs(min_velocity)) ){
                            lsf[i] = -lsf_top[i];
                        }else{
                            lsf[i] = lsf_top[i];
                        }
                    }
                    reset_lsf(100);
                    set_density();
                    set_volume();
                    volume_test = volume;
                }
        
            }

            while ( (fabs(volume_test-volume_target)/volume_in) > 0.005 && (counter < 100) && (perv_max < 1.0) && (perv_min > 0.0) ){
                //trecho não executado na primeira iteração
                ++counter;
                perv = 0.5*(perv_min + perv_max);
                for (int i=0; i<numNode; ++i){
                    if ( (lsf_top[i] < 0) && (normalVelocity[i] < 0) && (fabs(normalVelocity[i]) > perv*fabs(min_velocity)) ){                    
                        lsf[i] = -lsf_top[i];                    
                    }else{
                        lsf[i] = lsf_top[i];
                    }
                }
                reset_lsf(100);
                set_density();
                set_volume();
                volume_test = volume;

                if (volume_test > volume_target){                    
                    perv_max = perv;
                }else{                    
                    perv_min = perv;
                }
            }        
        }
    }//end if min_vel <> 0

    reset_lsf(100);

}   

void lset::set_regularizationMatrix_CSR(){
        
    double coef1, coef2;
    int mat_coord;

    coef1 = factorToRegularizeNormalVelocity/pow(elementLengthX,2);
    coef2 = factorToRegularizeNormalVelocity/pow(elementLengthY,2);
    
    //COO format assembly
	int NZ(0);
    for (int row = 0; row < (nely+1); ++row){
        for (int col = 0; col < (nelx+1); ++col){
            
            mat_coord = col*(nely+1)+row;
            
            coo_val[NZ] = 1+2*coef1+2*coef2;						
			coo_row[NZ] = mat_coord;
			coo_col[NZ] = mat_coord; 
			coo_idx[NZ] = NZ;
			++NZ;

            if (row>0){
                coo_val[NZ] = -coef2;						
			    coo_row[NZ] = mat_coord;
			    coo_col[NZ] = mat_coord-1; 
			    coo_idx[NZ] = NZ;
			    ++NZ;                
            }else{
                coo_val[NZ] = -coef2;						
			    coo_row[NZ] = mat_coord;
			    coo_col[NZ] = mat_coord; 
			    coo_idx[NZ] = NZ;
			    ++NZ;
            }
            
            if (row<nely){
                coo_val[NZ] = -coef2;						
			    coo_row[NZ] = mat_coord;
			    coo_col[NZ] = mat_coord+1; 
			    coo_idx[NZ] = NZ;
			    ++NZ;                
            }else{
                coo_val[NZ] = -coef2;						
			    coo_row[NZ] = mat_coord;
			    coo_col[NZ] = mat_coord; 
			    coo_idx[NZ] = NZ;
			    ++NZ;
            }

            if (col>0){
                coo_val[NZ] = -coef1;						
			    coo_row[NZ] = mat_coord;
			    coo_col[NZ] = mat_coord-(nely+1); 
			    coo_idx[NZ] = NZ;
			    ++NZ;
            }else{
                coo_val[NZ] = -coef1;						
			    coo_row[NZ] = mat_coord;
			    coo_col[NZ] = mat_coord; 
			    coo_idx[NZ] = NZ;
			    ++NZ;
            }

            if (col<nelx){
                coo_val[NZ] = -coef1;						
			    coo_row[NZ] = mat_coord;
			    coo_col[NZ] = mat_coord+(nely+1); 
			    coo_idx[NZ] = NZ;
			    ++NZ;                
            }else{
                coo_val[NZ] = -coef1;						
			    coo_row[NZ] = mat_coord;
			    coo_col[NZ] = mat_coord; 
			    coo_idx[NZ] = NZ;
			    ++NZ;                
            }
        }
    }
    
    //Assembly CSR format
    class sort_indices{
   		private:
    		int* mparr;
   		public:
     		sort_indices(int* parr) : mparr(parr) {}
     		bool operator()(int i, int j) const { return mparr[i]<mparr[j]; }
	};

	std::sort(coo_idx, coo_idx+NZ, sort_indices(coo_row));

    for (int i=0; i<NZ; ++i){
		csr_val[i] = coo_val[coo_idx[i]];
		csr_col[i] = coo_row[coo_idx[i]];
	}
    
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
       
	for (int i=1; i<counter; ++i){
		std::sort(coo_idx+csr_row[i-1], coo_idx+csr_row[i], sort_indices(coo_col));		
	}
	
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

void lset::solve_regularizationSystem_PCG(){
    /*---------------------------------------------------------------------------------------------
	Linear system solver using Left-preconditioned Conjugated Gradient
	Normal Residue (CGNR).	
	---------------------------------------------------------------------------------------------*/
	//matrix-vector multiplication
	int k1;
	int k2;
	int nummaxit(15000);//100//1000
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
	double stopCriteria(1e-15);//1e-6//1e-10//1e-9
	bool convergenceOK;

	convergenceOK = false;
	
	// ---windows only #pragma omp parallel for private(k1,k2,prod1) reduction(+:norm_b)
	for (int i=0; i<csr_rowsize; ++i){
		k1 = csr_row[i];
		k2 = csr_row[i+1];
		prod1 = 0;		
		for (int j=k1; j<k2;++j){
		    prod1 += csr_val[j]*normalVelocity_reg[csr_col[j]];
		}
		w[i] = prod1;
		norm_b += RHS_normalVelocity[i]*RHS_normalVelocity[i];
		r[i] = RHS_normalVelocity[i] - prod1;
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
	    normalVelocity_reg[i] += alp*p[i];
	    r[i] -= alp*w[i];
	    z[i] = cond_J[i]*r[i];
	    prod_zr2 += z[i]*r[i];
	    norm_u += normalVelocity_reg[i]*normalVelocity_reg[i];		
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
			normalVelocity_reg[i] += alp*p[i];
			r[i] -= alp*w[i];
			z[i] = cond_J[i]*r[i];
			prod_zr2 += z[i]*r[i];
			norm_u += normalVelocity_reg[i]*normalVelocity_reg[i];
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
		std::cout << " - lset PCG Convergence ERROR";
        solverNotConverged = true;
	}	
}

void lset::set_normalVelocity(){
    double aux_dv;
    aux_dv = elementLengthX*elementLengthY;

    for (int i=0; i<numNode; ++i){
        normalVelocity[i] = nodeCompliance[i]/aux_dv - lagrangeMultiplier_V - (1/penalizationParameter_V)*(volume-volumeTarget);       
    }

    regularize_normalVelocity();
    
}

void lset::regularize_normalVelocity(){

    double RegularizationMatrix_coef1, RegularizationMatrix_coef2;
    double aux_factorToAvoidZero;
    double dsx, dsy;

    int counter,mat_coord;

    aux_factorToAvoidZero = pow(factorToAvoidZero,2);

    RegularizationMatrix_coef1 = factorToRegularizeNormalVelocity/pow(elementLengthX,2);
    RegularizationMatrix_coef2 = factorToRegularizeNormalVelocity/pow(elementLengthY,2);
         
    for (int i=0; i<numNode; ++i){
        lsf_mirror[i] = lsf[i]/sqrt(pow(lsf[i],2)+aux_factorToAvoidZero);        
    }

    counter = 0;
    for (int col=0; col<(nelx+1); ++col){
        for (int row=0; row<(nely+1); ++row){

            dsx = 0.5*(set_dxp_neumannGrad(row,col)+set_dxm_neumannGrad(row,col));
            dsy = 0.5*(set_dyp_neumannGrad(row,col)+set_dym_neumannGrad(row,col));
                        
            RHS_normalVelocity[counter] = -normalVelocity[row*(nelx+1)+col]*0.5*sqrt(pow(dsx,2)+pow(dsy,2));
            ++counter;
        
        }
    }//end of main loop

    solve_regularizationSystem_PCG();
    counter = 0;
    for (int col=0; col<(nelx+1); ++col){
        for (int row=0; row<(nely+1); ++row){
            normalVelocity[row*(nelx+1)+col] = -normalVelocity_reg[counter];
            ++counter;
        }
    }    

}

void lset::set_normalVelocity_top(){

    for (int i=0; i<numNode; ++i){
        normalVelocity[i] = nodeTopGrad[i] - M_PI*(lagrangeMultiplier_V + (1/penalizationParameter_V)*(volume-volumeTarget));        
    }

}

void lset::removeLsfZeros(){

    srand( (unsigned)time( NULL ) );    

    for (int i=0; i<numNode; ++i){
        if (lsf[i]==0){
            lsf[i] = 1e-6*(2*round(static_cast <double> (rand()) / static_cast <double> (RAND_MAX)) -1);            
        }
    }

}

double lset::linearInterp(){

    double interpolated_density = 0;

    double sum_triangle_lsf;
    int sum_triangle_lsf_sign;
    
    elementLsf[4] = elementLsf[0]+elementLsf[1]+elementLsf[2]+elementLsf[3];
                
    if (elementLsf[4] == 0){
        elementLsf[4] = 1e-6*(2*round(static_cast <double> (rand()) / static_cast <double> (RAND_MAX)) -1);        
    }
    
    elementTriangleLsf[0] = elementLsf[4];
    elementTriangleLsf[1] = elementLsf[0];
    elementTriangleLsf[2] = elementLsf[1];
    elementTriangleLsf[3] = elementLsf[4];
    elementTriangleLsf[4] = elementLsf[1];
    elementTriangleLsf[5] = elementLsf[2];
    elementTriangleLsf[6] = elementLsf[4];
    elementTriangleLsf[7] = elementLsf[2];
    elementTriangleLsf[8] = elementLsf[3];
    elementTriangleLsf[9] = elementLsf[4];
    elementTriangleLsf[10] = elementLsf[3];
    elementTriangleLsf[11] = elementLsf[0];                  

    for (int row = 0; row<4; ++row){

        sum_triangle_lsf = 0;
        sum_triangle_lsf_sign = 0;                    
        
        for (int col=0; col<3; ++col){
            sum_triangle_lsf += elementTriangleLsf[row*3+col];
            sum_triangle_lsf_sign += sign(elementTriangleLsf[row*3+col]);            
        }

        if (sum_triangle_lsf_sign == -3){

            interpolated_density += 0.25;
            
                    
        }else if(sum_triangle_lsf_sign == 3){
                    
            interpolated_density += 0.25*ersatzMaterial;
            
                    
        }else{

            double dot_AB, dot_BB, dot_AC, dot_CC;
            double f1, f2;

            dot_AB = 0;
            dot_BB = 0;
            dot_AC = 0;
            dot_CC = 0;
            
            for (int col=0; col<3; ++col){

                lsf_low_aux[col] = 0;
                lsf_mid_aux[col] = 0;
                lsf_high_aux[col] = 0;

                if (sign(elementTriangleLsf[row*3+col]) != sign(sum_triangle_lsf)){ 
                    
                    switch(col){
                        case 0:
                            lsf_low_aux[col] = elementTriangleLsf[row*3+2];
                            lsf_mid_aux[col] = elementTriangleLsf[row*3];
                            lsf_high_aux[col] = elementTriangleLsf[row*3+1];
                            break;
                        case 1:
                            lsf_low_aux[col] = elementTriangleLsf[row*3];
                            lsf_mid_aux[col] = elementTriangleLsf[row*3+1];
                            lsf_high_aux[col] = elementTriangleLsf[row*3+2];
                            break;

                        case 2:
                    
                            lsf_low_aux[col] = elementTriangleLsf[row*3+1];
                            lsf_mid_aux[col] = elementTriangleLsf[row*3+2];
                            lsf_high_aux[col] = elementTriangleLsf[row*3];
                            break;
                    }
                                                          
                }
                
                dot_AB += lsf_mid_aux[col]*(lsf_mid_aux[col] - lsf_high_aux[col]);
                dot_BB += pow(lsf_mid_aux[col] - lsf_high_aux[col],2);
                dot_AC += lsf_mid_aux[col]*(lsf_mid_aux[col] - lsf_low_aux[col]);
                dot_CC += pow(lsf_mid_aux[col] - lsf_low_aux[col],2);

            }            

            f1 = dot_AB/dot_BB;
            f2 = dot_AC/dot_CC;       

            if (sum_triangle_lsf_sign == 1){
                interpolated_density += 0.25*(1-f1*f2)*ersatzMaterial+0.25*f1*f2;
                                 
            }else{
                interpolated_density += 0.25*(1-f1*f2)+0.25*f1*f2*ersatzMaterial;
                                
            }  
               
        }        

    }

    return interpolated_density;

}

double lset::set_dxm(int row, int col){    
    double dxm;

    if (col != 0){
        dxm = (lsf_mirror[row*(nelx+1)+col]-lsf_mirror[row*(nelx+1)+(col-1)])/elementLengthX;
    }else{
        dxm = 0;
    }    

    return dxm;
}

double lset::set_dxm_neumannGrad(int row, int col){    
    double dxm_ng;

    if (col != 0){
        dxm_ng = (lsf_mirror[row*(nelx+1)+col]-lsf_mirror[row*(nelx+1)+(col-1)])/elementLengthX;
    }else{
        dxm_ng = (lsf_mirror[row*(nelx+1)+(col+1)]-lsf_mirror[row*(nelx+1)+col])/elementLengthX;
    }    

    return dxm_ng;
}

double lset::set_dxmxm(int row, int col){
    double dxmxm;

    if (col != 0 && col!= 1){
        dxmxm = (lsf_mirror[row*(nelx+1)+col]-2*lsf_mirror[row*(nelx+1)+(col-1)]+lsf_mirror[row*(nelx+1)+(col-2)])/pow(elementLengthX,2);
    }else if(col == 0){
        dxmxm = 0;
    }else{
        dxmxm = (lsf_mirror[row*(nelx+1)+col]-lsf_mirror[row*(nelx+1)+(col-1)])/pow(elementLengthX,2);
    }   

    return dxmxm;
}

double lset::set_dxmxm_neumannGrad(int row, int col){
    double dxmxm_ng;

    if (col != 0 && col!= 1){
        dxmxm_ng = (lsf_mirror[row*(nelx+1)+col]-2*lsf_mirror[row*(nelx+1)+(col-1)]+lsf_mirror[row*(nelx+1)+(col-2)])/pow(elementLengthX,2);
    }else if(col == 0){
        dxmxm_ng = 0;
    }else{
        dxmxm_ng = 0;
    }   

    return dxmxm_ng;
}

double lset::set_dxp(int row, int col){
    double dxp;

    if (col != nelx){
        dxp = (lsf_mirror[row*(nelx+1)+(col+1)]-lsf_mirror[row*(nelx+1)+col])/elementLengthX;
    } else{
        dxp = 0;
    }

    return dxp;
}

double lset::set_dxp_neumannGrad(int row, int col){
    double dxp_ng;

    if (col != nelx){
        dxp_ng = (lsf_mirror[row*(nelx+1)+(col+1)]-lsf_mirror[row*(nelx+1)+col])/elementLengthX;
    } else{
        dxp_ng = (lsf_mirror[row*(nelx+1)+col]-lsf_mirror[row*(nelx+1)+(col-1)])/elementLengthX;
    }

    return dxp_ng;
}

double lset::set_dxpxp(int row, int col){
    double dxpxp;

    if (col != nelx-1 && col != nelx){
        dxpxp = (lsf_mirror[row*(nelx+1)+(col+2)]-2*lsf_mirror[row*(nelx+1)+(col+1)]+lsf_mirror[row*(nelx+1)+col])/pow(elementLengthX,2);
    }else if (col == nelx-1){
        dxpxp = (lsf_mirror[row*(nelx+1)+col]-lsf_mirror[row*(nelx+1)+(col+1)])/pow(elementLengthX,2);
    }else{
        dxpxp = 0;
    }

    return dxpxp;
}

double lset::set_dxpxp_neumannGrad(int row, int col){
    double dxpxp_ng;

    if (col != nelx-1 && col != nelx){
        dxpxp_ng = (lsf_mirror[row*(nelx+1)+(col+2)]-2*lsf_mirror[row*(nelx+1)+(col+1)]+lsf_mirror[row*(nelx+1)+col])/pow(elementLengthX,2);
    }else if (col == nelx-1){
        dxpxp_ng = 0;
    }else{
        dxpxp_ng = 0;
    }

    return dxpxp_ng;
}

double lset::set_dxmxp(int row, int col){
    double dxmxp;

    if (col != 0 && col != nelx){
        dxmxp = (lsf_mirror[row*(nelx+1)+(col-1)]-2*lsf_mirror[row*(nelx+1)+col]+lsf_mirror[row*(nelx+1)+(col+1)])/pow(elementLengthX,2);
    }else if (col == 0){
        dxmxp = (lsf_mirror[row*(nelx+1)+(col+1)]-lsf_mirror[row*(nelx+1)+col])/pow(elementLengthX,2);
    }else{
        dxmxp = (lsf_mirror[row*(nelx+1)+(col-1)]-lsf_mirror[row*(nelx+1)+col])/pow(elementLengthX,2);
    }
    
    return dxmxp;
}

double lset::set_dxmxp_neumannGrad(int row, int col){
    double dxmxp_ng;

    if (col != 0 && col != nelx){
        dxmxp_ng = (lsf_mirror[row*(nelx+1)+(col-1)]-2*lsf_mirror[row*(nelx+1)+col]+lsf_mirror[row*(nelx+1)+(col+1)])/pow(elementLengthX,2);
    }else if (col == 0){
        dxmxp_ng = 0;
    }else{
        dxmxp_ng = 0;
    }
    
    return dxmxp_ng;
}

double lset::set_dym(int row, int col){
    
    double dym;

    if (row != 0){
        dym = (lsf_mirror[row*(nelx+1)+col]-lsf_mirror[(row-1)*(nelx+1)+col])/elementLengthY;
    }else{
        dym = 0;
    }    

    return dym;
}

double lset::set_dym_neumannGrad(int row, int col){

    double dym_ng;

    if (row != 0){
        dym_ng = (lsf_mirror[row*(nelx+1)+col]-lsf_mirror[(row-1)*(nelx+1)+col])/elementLengthY;
    }else{
        dym_ng = (lsf_mirror[(row+1)*(nelx+1)+col]-lsf_mirror[row*(nelx+1)+col])/elementLengthY;
    }       

    return dym_ng;
}

double lset::set_dymym(int row, int col){
    double dymym;

    if (row != 0 && row!= 1){
        dymym = (lsf_mirror[row*(nelx+1)+col]-2*lsf_mirror[(row-1)*(nelx+1)+col]+lsf_mirror[(row-2)*(nelx+1)+col])/pow(elementLengthY,2);
    }else if(row == 0){
        dymym = 0;
    }else{
        dymym = (lsf_mirror[row*(nelx+1)+col]-lsf_mirror[(row-1)*(nelx+1)+col])/pow(elementLengthY,2);
    }   

    return dymym;
}

double lset::set_dymym_neumannGrad(int row, int col){
    double dymym_ng;

    if (row != 0 && row!= 1){
        dymym_ng = (lsf_mirror[row*(nelx+1)+col]-2*lsf_mirror[(row-1)*(nelx+1)+col]+lsf_mirror[(row-2)*(nelx+1)+col])/pow(elementLengthY,2);
    }else if(row == 0){
        dymym_ng = 0;
    }else{
        dymym_ng = 0;
    }   

    return dymym_ng;
}

double lset::set_dyp(int row, int col){
    double dyp;

    if (row != nely){
        dyp = (lsf_mirror[(row+1)*(nelx+1)+col]-lsf_mirror[row*(nelx+1)+col])/elementLengthY;
    } else{
        dyp = 0;
    }

    return dyp;
}

double lset::set_dyp_neumannGrad(int row, int col){
    double dyp_ng;

    if (row != nely){
        dyp_ng = (lsf_mirror[(row+1)*(nelx+1)+col]-lsf_mirror[row*(nelx+1)+col])/elementLengthY;
    } else{
        dyp_ng = (lsf_mirror[(row)*(nelx+1)+col]-lsf_mirror[(row-1)*(nelx+1)+col])/elementLengthY;
    }

    return dyp_ng;
}

double lset::set_dypyp(int row, int col){
    double dypyp;

    if (row != nely-1 && row != nely){
        dypyp = (lsf_mirror[(row+2)*(nelx+1)+col]-2*lsf_mirror[(row+1)*(nelx+1)+col]+lsf_mirror[row*(nelx+1)+col])/pow(elementLengthY,2);
    }else if (row == nely-1){
        dypyp = (lsf_mirror[row*(nelx+1)+col]-lsf_mirror[(row+1)*(nelx+1)+col])/pow(elementLengthY,2);
    }else{
        dypyp = 0;
    }
    
    return dypyp;
}

double lset::set_dypyp_neumannGrad(int row, int col){
    double dypyp_ng;

    if (row != nely-1 && row != nely){
        dypyp_ng = (lsf_mirror[(row+2)*(nelx+1)+col]-2*lsf_mirror[(row+1)*(nelx+1)+col]+lsf_mirror[row*(nelx+1)+col])/pow(elementLengthY,2);
    }else if (row == nely-1){
        dypyp_ng = 0;
    }else{
        dypyp_ng = 0;
    }

    return dypyp_ng;
}

double lset::set_dymyp(int row, int col){
    double dymyp;

    if (row != 0 && row != nely){
        dymyp = (lsf_mirror[(row+1)*(nelx+1)+col]-2*lsf_mirror[row*(nelx+1)+col]+lsf_mirror[(row-1)*(nelx+1)+col])/pow(elementLengthY,2);
    }else if (row == 0){
        dymyp = (lsf_mirror[(row+1)*(nelx+1)+col]-lsf_mirror[row*(nelx+1)+col])/pow(elementLengthY,2);
    }else{
        dymyp = (lsf_mirror[(row-1)*(nelx+1)+col]-lsf_mirror[row*(nelx+1)+col])/pow(elementLengthY,2);
    }
    
    return dymyp;
}

double lset::set_dymyp_neumannGrad(int row, int col){
    double dymyp_ng;

    if (row != 0 && row != nely){
        dymyp_ng = (lsf_mirror[(row+1)*(nelx+1)+col]-2*lsf_mirror[row*(nelx+1)+col]+lsf_mirror[(row-1)*(nelx+1)+col])/pow(elementLengthY,2);
    }else if (row == 0){
        dymyp_ng = 0;
    }else{
        dymyp_ng = 0;
    }
    
    return dymyp_ng;
}

double lset::minmod(double arg1, double arg2){
    
    return static_cast <double>(std::max(0,sign(arg1)*sign(arg2))*sign(arg1))*std::min(std::abs(arg1),std::abs(arg2)); 
}

double lset::flux(double u1,double u2, double v1, double v2){
    
    return sqrt(pow(std::max(u1,0.0),2) + pow(std::min(u2,0.0),2) + pow(std::max(v1,0.0),2) + pow(std::min(v2,0.0),2));
}

void lset::set_timeStep(){

    double maxNormalVelocity;
		
	maxNormalVelocity = 0;
	for (int i = 0; i < numNode; ++i){
        if (std::abs(normalVelocity[i]) > maxNormalVelocity){
            maxNormalVelocity = std::abs(normalVelocity[i]);
        }
    }
    
    timeStep = CFL*0.5*std::min(elementLengthX,elementLengthY)/maxNormalVelocity;
}

void lset::set_curvature(){
    
    double gradx_lsf, grady_lsf, mag;
    double div_normalVector_x, div_normalVector_y;
    
    // ---windows only #pragma omp parallel for
    for (int i=0; i<numNode; ++i){
        lsf_mirror[i] = lsf[i];
    }
    
    // ---windows only #pragma omp parallel for private(gradx_lsf,grady_lsf,mag)
    for (int row=0; row<(nely+1); ++row){
        for (int col=0; col<(nelx+1); ++col){
            gradx_lsf = 0.5*(set_dxp_neumannGrad(row,col)+set_dxm_neumannGrad(row,col));
            grady_lsf = 0.5*(set_dyp_neumannGrad(row,col)+set_dym_neumannGrad(row,col));
            mag = sqrt(pow(gradx_lsf,2)+pow(grady_lsf,2)+pow(factorToAvoidZero,2));

            normalVector_x[row*(nelx+1)+col] = gradx_lsf/mag;
            normalVector_y[row*(nelx+1)+col] = grady_lsf/mag;
        }
    }
    
    // ---windows only #pragma omp parallel for
    for (int i=0; i<numNode; ++i){
        lsf_mirror[i] = normalVector_x[i];
        curvature[i] = 0;
    }

    // ---windows only #pragma omp parallel for private(div_normalVector_x)
    for (int row=0; row<(nely+1); ++row){
        for (int col=0; col<(nelx+1); ++col){
            div_normalVector_x = 0.5*(set_dxp_neumannGrad(row,col)+set_dxm_neumannGrad(row,col));
            curvature[row*(nelx+1)+col] += div_normalVector_x;
        }
    }
    
    // ---windows only #pragma omp parallel for
    for (int i=0; i<numNode; ++i){
        lsf_mirror[i] = normalVector_y[i];       
    }

    // ---windows only #pragma omp parallel for private(div_normalVector_y)
    for (int row=0; row<(nely+1); ++row){
        for (int col=0; col<(nelx+1); ++col){
            div_normalVector_y = 0.5*(set_dyp_neumannGrad(row,col)+set_dym_neumannGrad(row,col));
            curvature[row*(nelx+1)+col] += div_normalVector_y;
        }
    }
        
}

void lset::set_density(){

    int sum_node_lsf_sign;
                
    srand( (unsigned)time( NULL ) );

    removeLsfZeros();
    
    for (int i=0; i<(numElement); ++i){

        density[i] = 0;
        sum_node_lsf_sign = 0;

        for (int j=0; j<4; ++j){
            
            elementLsf[j] = lsf[element[i*8+(4+j)]-1];
            sum_node_lsf_sign += sign(elementLsf[j]);
            
        }       

        if (sum_node_lsf_sign == 4){

            density[i] = ersatzMaterial;

        }else if(sum_node_lsf_sign == -4){

            density[i] = 1.0;
            
        }else{         

            density[i] = linearInterp();     

        }
        
    }

}

void lset::set_volume(){
    volume = 0;
    for (int i=0; i<(numElement); ++i){
        volume += density[i];
    }
    volume = volume*elementLengthX*elementLengthY;    
}

void lset::set_volumeMax(){
    volumeMax = 0;
    for (int i=0; i<(numElement); ++i){
        ++volumeMax;
    }
    volumeMax = volumeMax*elementLengthX*elementLengthY;

    volumeTarget = volFrac_target*volumeMax;
}

void lset::set_volumeFraction(){
    volumeFraction = volume/volumeMax;
}

void lset::set_perimeter(){
    double factor_perimeter;
    
    double dsy,dsx;

    factor_perimeter = pow(factorToAvoidZero,2);
    
    perimeter = 0;

    // ---windows only #pragma omp parallel for
    for (int i=0; i<numNode; ++i){
        lsf_mirror[i] = lsf[i]/sqrt(pow(lsf[i],2)+factor_perimeter);
    }

    // ---windows only #pragma omp parallel for private(dsx,dsy) reduction(+:perimeter)
    for (int row=0; row<(nely+1); ++row){
        for (int col=0; col<(nelx+1); ++col){

            dsx = 0.5*(set_dxp_neumannGrad(row,col)+set_dxm_neumannGrad(row,col));
            dsy = 0.5*(set_dyp_neumannGrad(row,col)+set_dym_neumannGrad(row,col));
            perimeter += sqrt(pow(dsx,2)+pow(dsy,2));

        }
    }//end of main loop
    
    perimeter = 0.5*elementLengthX*elementLengthY*perimeter;

}

double lset::set_objectiveFunction(){

    set_volume();
    set_volumeFraction();
    set_perimeter();    
    
    double hi;
    hi = volume-volumeTarget;

    return compliance[0] + lagrangeMultiplier_V*hi + 0.5*(1/penalizationParameter_V)*pow(hi,2) + lagrangeMultiplier_P*perimeter;
    
}

void lset::set_seedHole(int iterId){
    
    if ((iterId+1)%solve_lsf_top_freq != 0){
        seedHole = false;
    }else{
        seedHole = true;
    }

}

void lset::set_lagrangeMultiplier_V(){
    
    penalizationParameter_V = std::max(0.9*penalizationParameter_V,0.1); 
    lagrangeMultiplier_V += (1/penalizationParameter_V)*(volume - volumeTarget);
    
}

// bool lset::check_topologyConvergence(int numIter, int lset_solutionAttempt){
bool lset::check_topologyConvergence(int numIter, int lset_solutionAttempt, int objectiveFunction_step){

    bool topologyConvergence;

    // <new>
    double tolerance(1e-3);
    double optimizationConvergenceCriteria;
    // <new>

    if (numIter >= num_iter_max){
        topologyConvergence = true;
    }else{        
        // if ( (fabs(volume - volumeTarget)/volumeMax < 0.01) && lset_solutionAttempt == 6){
        //     topologyConvergence = true;
        // }else{
        //     topologyConvergence = false;
        // }

        // <new>
        if (flag_initialTopology){
            if(objectiveFunction_step>2){
                optimizationConvergenceCriteria = (objectiveFunction[objectiveFunction_step] - objectiveFunction[objectiveFunction_step-1])/(1+objectiveFunction[objectiveFunction_step-1]);
                if ((fabs(volume - volumeTarget)/volumeMax < 0.01) && optimizationConvergenceCriteria < tolerance){
                    topologyConvergence = true;
                }else{
                    topologyConvergence = false;
                }
            }else{
                topologyConvergence = false;
            }
        }else{
            if ( (fabs(volume - volumeTarget)/volumeMax < 0.01) && lset_solutionAttempt == 6){
                topologyConvergence = true;
            }else{
                topologyConvergence = false;
            }
        }
        // <new>

    }

    return topologyConvergence; 

}

void lset::log_volumeFraction(){

    std::ofstream myvolumefractionfile;

    myvolumefractionfile.open (dirPath + "\\volumeFraction.txt",std::ios_base::app);

    myvolumefractionfile << std::setprecision(5) << volumeFraction << " ";

    myvolumefractionfile.close();

}

void lset::log_objectiveFunction(double objectiveFunction){

    std::ofstream myobjectivefunctionfile;

    myobjectivefunctionfile.open (dirPath + "\\objectiveFunction.txt",std::ios_base::app);

    myobjectivefunctionfile << std::setprecision(15) << objectiveFunction << " ";

    myobjectivefunctionfile.close();

}

void lset::log_iterationNumber(int iterNumber){

    std::ofstream myiterationnumberfile;

    myiterationnumberfile.open (dirPath + "\\iterationCounting.txt",std::ios_base::app);

    myiterationnumberfile << iterNumber << " ";

    myiterationnumberfile.close();

}

void lset::log_compliance(){

    std::ofstream mycompliancefile;

    mycompliancefile.open (dirPath + "\\compliance.txt",std::ios_base::app);

    mycompliancefile << std::setprecision(15) << compliance[0] << " ";

    mycompliancefile.close();

}

void lset::log_density_hist(){
     
    std::ofstream mydensityfile;
    
    mydensityfile.open (dirPath + "\\density-hist.txt",std::ios_base::app);

    for (int i=0; i<(numElement); ++i){
        mydensityfile << std::setprecision(15) << density[i] << " ";        
    } 
    mydensityfile << "\n";
    mydensityfile.close();
}

void lset::log_lsf_hist(){
     
    std::ofstream mylsffile;
    
    mylsffile.open (dirPath + "\\lsf-hist.txt",std::ios_base::app);

    for (int i=0; i<(numNode); ++i){
        mylsffile << std::setprecision(15) << lsf[i] << " ";        
    } 
    mylsffile << "\n";
    mylsffile.close();
}

void lset::log_data_finalresult(int iterNumber, bool convergenceOK, double simuTime){

    log_simulationData(iterNumber, convergenceOK,simuTime);
    //log_density();
}

void lset::log_simulationData(int iterNumber, bool convergenceOk, double simuTime){

    std::ofstream mysimudatafile;
    mysimudatafile.open (dirPath + "\\_simudata.txt");

    mysimudatafile << meshFilePath << std::endl;

    //simulation data
    mysimudatafile << "convergence acchieved: " << convergenceOk << std::endl;
    mysimudatafile << "number of iterations: " << iterNumber << std::endl;
    mysimudatafile << "Time elapsed during simulation (sec): " << simuTime << std::endl;    

    //mesh data
    mysimudatafile << "number of elements x-direction: " << nelx << std::endl;
    mysimudatafile << "number of elements y-direction: " << nely << std::endl;
    mysimudatafile << "element lenght x-direction: " << elementLengthX << std::endl;
    mysimudatafile << "element lenght y-direction: " << elementLengthY << std::endl;
    mysimudatafile << "design domain lenght x-direction: " << lengthX << std::endl;
    mysimudatafile << "design domain lenght y-direction: " << lengthY << std::endl;
    
    //level set data
    mysimudatafile << "number of holes x-direction: " << hole_numx << std::endl;    
    mysimudatafile << "number of holes y-direction: " << hole_numy << std::endl;
    mysimudatafile << "hole radius: " << hole_r << std::endl;
    mysimudatafile << "number of steps to reset lsf: " << num_step_lsf_reset << std::endl;
    mysimudatafile << "frequence to reset lsf: " << reset_lsf_freq << std::endl;
    mysimudatafile << "initial number of steps to solve lsf: " << num_step_lsf_solve << std::endl;
    mysimudatafile << "max iteration steps number: " << num_iter_max << std::endl;
    mysimudatafile << "lagrange multiplier volume: " << lagrangeMultiplier_V << std::endl;
    mysimudatafile << "penalization parameter volume: " << penalizationParameter_V << std::endl;
    mysimudatafile << "lagrange multiplier perimeter: " << lagrangeMultiplier_P << std::endl;
    mysimudatafile << "allowance beetwing iteration - advection: " << allowance_adv << std::endl;
    mysimudatafile << "allowance beetwing iteration - topology: " << allowance_top << std::endl;
    mysimudatafile << "frequence to topology gradient calculation: " << solve_lsf_top_freq << std::endl;
    mysimudatafile << "fraction to remove on each top grad computation: " << frac_to_remove << std::endl;
    mysimudatafile << "target volume fraction: " << volFrac_target << std::endl;
      

    mysimudatafile.close();

}

void lset::remove_files(std::string filename){
    
    std::string density_name = "repo/cap4/" + filename + "_density.csv";
    std::string lsf_name = "repo/cap4/" + filename + "_lsf.csv";
    std::string fobcurve_name = "repo/cap4/" + filename + "_fobcurve.csv";

    const char *str2char_density_name = density_name.c_str();
    const char *str2char_lsf_name = lsf_name.c_str();
    const char *str2char_fobcurve_name = fobcurve_name.c_str();

    std::remove(str2char_density_name);
    std::remove(str2char_lsf_name);
    std::remove(str2char_fobcurve_name);
    
}

void lset::log2files(std::string filename, int iter, double fob, double iter_compliance, double volfrac){

    std::string density_name = "repo/cap4/" + filename + "_density.csv";
    std::string lsf_name = "repo/cap4/" + filename + "_lsf.csv";
    std::string fobcurve_name = "repo/cap4/" + filename + "_fobcurve.csv";

    std::ofstream mydensityfile;
    mydensityfile.open (density_name,std::ios_base::app);
    for (int i=0; i<(numElement-1); ++i){
        mydensityfile << std::setprecision(17) << density[i] << ",";        
    } 
    mydensityfile << density[numElement-1] << "\n";
    mydensityfile.close();

    std::ofstream mylsffile;
    mylsffile.open (lsf_name,std::ios_base::app);
    for (int i=0; i<(numNode-1); ++i){
        mylsffile << std::setprecision(17) << lsf[i] << ",";        
    } 
    mylsffile << lsf[numNode-1] << "\n";
    mylsffile.close();

    std::ofstream myfobcurvefile;
    mylsffile.open (fobcurve_name,std::ios_base::app);
    mylsffile << iter << ",";
    mylsffile << std::setprecision(17) << fob << ","; 
    mylsffile << std::setprecision(17) << iter_compliance << ","; 
    mylsffile << std::setprecision(17) << volfrac << "\n"; 
    mylsffile.close(); 

}

void lset::log2files_red(std::string filename, int iter, double fob, double iter_compliance, double volfrac){

    std::string density_name_red = "repo/cap4/" + filename + "_density_red.csv";
    std::string lsf_name_red = "repo/cap4/" + filename + "_lsf_red.csv";
    std::string fobcurve_name_red = "repo/cap4/" + filename + "_fobcurve_red.csv";

    std::ofstream mydensityfile;
    mydensityfile.open (density_name_red,std::ios_base::app);
    for (int i=0; i<(numElement-1); ++i){
        mydensityfile << std::setprecision(17) << density[i] << ",";        
    } 
    mydensityfile << density[numElement-1] << "\n";
    mydensityfile.close();

    std::ofstream mylsffile;
    mylsffile.open (lsf_name_red,std::ios_base::app);
    for (int i=0; i<(numNode-1); ++i){
        mylsffile << std::setprecision(17) << lsf[i] << ",";        
    } 
    mylsffile << lsf[numNode-1] << "\n";
    mylsffile.close();

    std::ofstream myfobcurvefile;
    mylsffile.open (fobcurve_name_red,std::ios_base::app);
    mylsffile << iter << ",";
    mylsffile << std::setprecision(17) << fob << ","; 
    mylsffile << std::setprecision(17) << iter_compliance << ","; 
    mylsffile << std::setprecision(17) << volfrac << "\n"; 
    mylsffile.close(); 

}

void lset::log_density(std::string filename){
    
    std::string density_name = "repo/cap4/" + filename + "_density.csv";

    std::ofstream mydensityfile;
    
    //mydensityfile.open ("result-repo/density.csv");
    //mydensityfile.open ("repo/cap4/bridge_density.csv",std::ios_base::app);
    mydensityfile.open (density_name,std::ios_base::app);

    for (int i=0; i<(numElement-1); ++i){
        mydensityfile << std::setprecision(17) << density[i] << ",";        
    } 
    mydensityfile << density[numElement-1] << "\n";
    mydensityfile.close();
}

void lset::log_lsf(std::string filename){
    
    std::string lsf_name = "repo/cap4/" + filename + "_density.csv";

    std::ofstream mylsffile;
    
    //mydensityfile.open ("result-repo/density.csv");
    //mylsffile.open ("repo/cap4/bridge_lsf.csv",std::ios_base::app);
    mylsffile.open (lsf_name,std::ios_base::app);


    for (int i=0; i<(numNode-1); ++i){
        mylsffile << std::setprecision(17) << lsf[i] << ",";        
    } 
    mylsffile << lsf[numNode-1] << "\n";
    mylsffile.close();
}

//=================================================================================================
// Public member functions
//=================================================================================================

void lset::optimize_topology(fem &femmodel){

    int numIter(0);
    int lset_solutionAttempt(1);
    int auxToFindMax;
    double objectiveFunction_test;
    double allowance;
    int iterNumber(0);
    // <new>
    int objectiveFunction_step(0); 
    // <new>
    double simuTime(0);    
    
    bool topologyConvergence;
    bool convergenceOK;

    topologyConvergence = false;
    convergenceOK = false;
    CFL = 1.0;

    femmodel.set_femAnalysis();
    solverNotConverged = femmodel.get_solverNotConverged();
    // <new>
    for (int i = 0; i<(numElement);++i){
        elementCompliance_init[i] = elementCompliance[i];
        elementTopGrad_init[i] = elementTopGrad[i];
    }
    // <new>
    
    if (flag_topologicalGradient == true){
        set_normalVelocity_top();                
    }else{
        set_normalVelocity();
    }    
    
    // objectiveFunction = set_objectiveFunction();
    // <new>
    objectiveFunction[objectiveFunction_step] = set_objectiveFunction();
    // <new>

    while (topologyConvergence == false && solverNotConverged == false){

        for (int i=0; i<numNode; ++i){
            lsf_tmp[i] = lsf[i]; 
        }

        if (flag_topologicalGradient == true){
    
            set_seedHole(numIter);
            if ((seedHole == true) || (numIter == 0)){
                solve_lsf_top();
            }else{
                set_timeStep();
                solve_lsf(num_step_lsf_solve);
            }
        
        }else{        
            set_timeStep();
            solve_lsf(num_step_lsf_solve);
        }
              
        set_density();
                
        femmodel.set_femAnalysis();
        solverNotConverged = femmodel.get_solverNotConverged();
                
        objectiveFunction_test = set_objectiveFunction();

        //std::cout << "step " << numIter+1 << " of " << num_iter_max;
        // std::cout << " | objFunction = " << objectiveFunction;
        // <new>
        //std::cout << " | objFunction = " << objectiveFunction[objectiveFunction_step];
        // <new>
        //std::cout << " | objFunction_test = " << objectiveFunction_test;
        //std::cout << " | HJi = " << num_step_lsf_solve;
        //std::cout << " | HJa = " << lset_solutionAttempt;
        //std::cout << " | Volume = " << volumeFraction << std::endl;        


        //set allowance
        if (flag_topologicalGradient == true){

            if ((seedHole == true) || (numIter == 0)){
                allowance = allowance_top;
            }else{
                allowance = allowance_adv;
            }
        
        }else{
            
            allowance = allowance_adv;            
        }

        // if (numIter+1 >= 0.75*num_iter_max ){
        //     allowance = 0;
        //     flag_topologicalGradient = false;
        // }
        // <new>
        if (flag_initialTopology){
            if (numIter+1 >= 0.10*num_iter_max){
                allowance = 0;
                flag_topologicalGradient = false;
            }
        }else{
            if (numIter+1 >= 0.75*num_iter_max){
                allowance = 0;
                flag_topologicalGradient = false;
            }
        }
        // <new>

        // if (objectiveFunction_test<=objectiveFunction*(1+allowance)){
        // <new>
        if (objectiveFunction_test<=objectiveFunction[objectiveFunction_step]*(1+allowance)){
        // <new>
            ++numIter;
            // objectiveFunction = objectiveFunction_test;
            // <new>
            ++objectiveFunction_step;
            objectiveFunction[objectiveFunction_step] = objectiveFunction_test;
            // <new>

            if (flag_topologicalGradient == true){                
                
                set_seedHole(numIter);
                if ((seedHole == true) || (numIter%solve_lsf_top_freq == 0) || (numIter == 1)){                    
                    femmodel.set_femAnalysis();
                    solverNotConverged = femmodel.get_solverNotConverged();
                    if (seedHole == true){
                        set_normalVelocity_top();
                    }else{
                        set_normalVelocity();                                                
                    }                
                
                }else{                    
                    set_normalVelocity();                                  
                }
            }else{                
                set_normalVelocity();            
            }
                        
            auxToFindMax = 1.1*num_step_lsf_solve;
            num_step_lsf_solve = std::min(10,std::max(auxToFindMax,num_step_lsf_solve+1));
            CFL = 1.0;
            lset_solutionAttempt = 1;

            //update lagrangian multipliers
            set_lagrangeMultiplier_V();

            ++iterNumber;
                    
        }else{//objective_test is higher then objective
            ++numIter;
            
            for (int i=0; i<numNode; ++i){
                lsf[i] = lsf_tmp[i];                
            }
            
            set_density();

            if (flag_topologicalGradient == true){
                
                set_seedHole(numIter);
                if ((seedHole == true) || (numIter%solve_lsf_top_freq == 0)){
                    femmodel.set_femAnalysis();
                    solverNotConverged = femmodel.get_solverNotConverged();
                    if (seedHole == true){
                        set_normalVelocity_top();
                    }else{
                        set_normalVelocity();                        
                    }
                }else{
                    num_step_lsf_solve = 0.5*num_step_lsf_solve;
                    if (num_step_lsf_solve==0){
                        num_step_lsf_solve = 1;
                        CFL = 0.5*CFL;
                    }
                    ++lset_solutionAttempt;                    
                }

                // topologyConvergence = check_topologyConvergence(numIter, lset_solutionAttempt);
                // <new>
                topologyConvergence = check_topologyConvergence(numIter, lset_solutionAttempt,objectiveFunction_step);
                // <new>

            }else{
                
                num_step_lsf_solve = 0.5*num_step_lsf_solve;
                if (num_step_lsf_solve==0){
                    num_step_lsf_solve = 1;
                    CFL = 0.5*CFL;
                }
                ++lset_solutionAttempt;                
                // topologyConvergence = check_topologyConvergence(numIter, lset_solutionAttempt);
                // <new>
                topologyConvergence = check_topologyConvergence(numIter, lset_solutionAttempt,objectiveFunction_step);
                // <new>
            }

        }     

    }//end while loop    
    
    if(numIter == num_iter_max || solverNotConverged == true){
        convergenceOK = false;
        // std::cout << " Optimization procedure reached maximum allowable interation before convergence" << std::endl;
        std::cout << " - xxx BAD xxx - convergence not acchieved";
        sim_converged = 0;
        // <new>
        num_iteration = numIter;
        // <new>
    }else{
        convergenceOK = true;
        std::cout << " - <<< OK >>> convergence acchieved";
        // std::cout << " Topology optimization completed" << std::endl;
        //log_density();
        sim_converged = 1;
        // <new>
        num_iteration = numIter;
        // <new>
    }  

}

void lset::get_initial_compliance_and_topgrad(fem &femmodel){
    
    femmodel.set_femAnalysis();
    // <data saving test>
    for (int i = 0; i<(numElement);++i){
        elementCompliance_init[i] = elementCompliance[i];
        elementTopGrad_init[i] = elementTopGrad[i];
    }

}
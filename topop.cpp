// code based on Allaire for topology optimization with level set
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include <cmath>
#include <chrono>
#include <ctime>
#include <iostream>
#include <string>
#include <time.h>

#include <chrono>

#include "batch.h"
#include "mesh.h"
#include "fem.h"
#include "lset.h"


//int main(){
int main(int argc, char** argv){
    
    optInput lsetInput;
    //<new>
    optMode lsetMode;
    //<new>
    int batch_size(0);
    std::string line;
    std::string simulationId;
    std::string meshFilePath;
    std::string loadFilePath;
    std::string boundaryFilePath;
    std::string *lset_params;
    int simulationFlag;
    int simulationConvergence;
    double simulationTime;
    double *final_density;
    double *initial_elementCompliance;
    double *initial_elementTopGrad;
    std::string final_density_name;
    std::string initial_compliance_name;
    std::string initial_topgrad_name;
    //<new>
    int num_iteration;
    double final_compliance;
    //<new>

    std::chrono::steady_clock::time_point begin; 
    std::chrono::steady_clock::time_point end; 

    //batch mybatch("repo/","batch_file_05.csv");
    batch mybatch("repo/",argv[1]);

    batch_size = mybatch.get_batch_size();

    std::cout << std::endl <<"Starting batch processing" << std::endl;

    for (int row = 0; row < batch_size; ++row){

        if (row == 0){
            std::cout << std::endl << row + 1 << " of " << batch_size;
        }else{
            std::cout << std::endl << row + 1 << " of " << batch_size;
        }
        
        simulationId = mybatch.get_simulationId(row);
        meshFilePath = mybatch.get_meshFilePath(row);
        loadFilePath = mybatch.get_loadFilePath(row);
        boundaryFilePath = mybatch.get_boundaryFilePath(row);
        lset_params = mybatch.get_lset_params(row);
        simulationFlag = mybatch.get_simulationFlag(row);

        lsetInput.hole_num_x = stoi(lset_params[0]);
        lsetInput.hole_num_y = stoi(lset_params[1]);
        lsetInput.hole_radius = stold(lset_params[2]);
        lsetInput.reset_lsf_freq = stoi(lset_params[3]);
        lsetInput.reset_lsf_numStep = stoi(lset_params[4]);
        lsetInput.solve_lsf_numStep_initial = stoi(lset_params[5]);
        lsetInput.topologyOptimization_numIter_max = stoi(lset_params[6]);
        lsetInput.lagrangeMultiplier_V = stoi(lset_params[7]);
        lsetInput.penalizationParameter_V = stoi(lset_params[8]);
        lsetInput.lagrangeMultiplier_P = stold(lset_params[9]);
        lsetInput.allowance_adv = stold(lset_params[10]);
        lsetInput.solve_lsf_top_freq = stoi(lset_params[11]);
        lsetInput.solve_lsf_top_fracToRemove = stold(lset_params[12]);
        lsetInput.allowance_top = stold(lset_params[13]);
        lsetInput.volume_fraction_target = stold(lset_params[14]);

        //<new>
        lsetMode.optMode = mybatch.get_lsetMode(row);
        lsetMode.inputDensityFile = mybatch.get_inputDensity(row);
        //<new>

        // lsetInput.hole_num_x = 4;
        // lsetInput.hole_num_y = 2;
        // lsetInput.hole_radius = 0.5;
        // lsetInput.reset_lsf_freq = 2;
        // lsetInput.reset_lsf_numStep = 2;
        // lsetInput.solve_lsf_numStep_initial = 30;
        // lsetInput.topologyOptimization_numIter_max = 300;
        // lsetInput.lagrangeMultiplier_V = 10;
        // lsetInput.penalizationParameter_V = 100;
        // lsetInput.lagrangeMultiplier_P = 0.1;
        // lsetInput.allowance_adv = 0.01;
        // lsetInput.solve_lsf_top_freq = 6;
        // lsetInput.solve_lsf_top_fracToRemove = 0.04;
        // lsetInput.allowance_top = 0.10;
        // lsetInput.volume_fraction_target = 0.50;

        
        if (simulationFlag == 0){
            
            mesh mymesh(meshFilePath, boundaryFilePath, loadFilePath);
            fem fea(mymesh);
            lset top(mymesh, fea, lsetInput, lsetMode);
            
            begin = std::chrono::steady_clock::now();
            top.optimize_topology(fea);
            end = std::chrono::steady_clock::now();
            
            simulationTime = std::chrono::duration_cast<std::chrono::seconds>(end - begin).count();
            
            simulationConvergence = top.get_sim_converged();

            //<new - compliance and topgrad names for files saving>
            initial_compliance_name = simulationId + "_compliance.csv";
            initial_topgrad_name = simulationId + "_topgrad.csv";

            initial_elementCompliance = top.get_elementCompliance_init();
            initial_elementTopGrad = top.get_elementTopGrad_init();
            
            mybatch.set_simulationInitialElementCompliance(initial_elementCompliance,mymesh.get_numElement(),initial_compliance_name);
            mybatch.set_simulationInitialElementTopGrad(initial_elementTopGrad,mymesh.get_numElement(),initial_topgrad_name);
            //<new>
        
            final_density_name = simulationId + "_top.csv";
            final_density = top.get_density();
            mybatch.set_simulationFinalTopologyFile(final_density, mymesh.get_numElement(), final_density_name);

            num_iteration = top.get_num_iteration();
            final_compliance = top.get_compliance();
            
            mybatch.update_batch_file(row, 
                                      simulationConvergence, 
                                      simulationTime,
                                      final_density_name,
                                      num_iteration,
                                      final_compliance);
            
 
        }

    }
    
    std::cout << std::endl <<"Batch processing completed" << std::endl;
    return 0;
}
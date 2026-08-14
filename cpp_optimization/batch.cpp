#include "batch.h"
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <cmath>

batch::batch(std::string batch_file_path, std::string batch_name):
    batch_root(batch_file_path),
    batch_path(batch_root+batch_name),
    mesh_root(batch_root+"meshes/"),
    bdry_root(batch_root+"boundaries/"),
    load_root(batch_root+"loads/"),
    params_root(batch_root+"params/"),
    final_topology_root(batch_root+"final_topologies/"),
    initial_compliance_root(batch_root+"initial_compliance/"),
    initial_topgrad_root(batch_root+"initial_topgrad/"),
    batch_size(0),
    lset_params(nullptr),
    batch_dataframe(nullptr){

        set_batch_size();
        lset_params = new std::string [15]();

        set_batch_dataframe();
    }

batch::~batch(){
    delete [] batch_dataframe;
}

//==================================================================================================
//private member functions
//==================================================================================================
void batch::set_batch_size(){
	
    std::ifstream input_batch_file;
    std::string line;

    input_batch_file.open(batch_path, std::ios::in);
    if (input_batch_file.good()){
        while (!input_batch_file.eof()){
            getline(input_batch_file, line);
            ++batch_size;
        }
        input_batch_file.close();
    }
    else{
        std::cout << "error: check batch_file" << std::endl;
    }
}

void batch::set_batch_dataframe(){
    
    std::ifstream input_batch_file;
    std::string s1;
    std::string s2;
    std::string s3;
    std::string s4;
    std::string s5;
    std::string s6;
    std::string s7;
    std::string s8;
    std::string s9;
    // <new>
    std::string s10;
    std::string s11;
    std::string s12;
    std::string s13;
    // <new>

    batch_dataframe = new batchData [batch_size]();
    
    input_batch_file.open(batch_path, std::ios::in);

    for (int row = 0 ; row < batch_size; ++row){
        std::getline(input_batch_file, s1,',');
        std::getline(input_batch_file, s2,',');
        std::getline(input_batch_file, s3,',');
        std::getline(input_batch_file, s4,',');
        std::getline(input_batch_file, s5,',');
        std::getline(input_batch_file, s6,',');
        std::getline(input_batch_file, s7,',');
        std::getline(input_batch_file, s8,',');
        //std::getline(input_batch_file, s9);
        //<new>
        std::getline(input_batch_file, s9,',');
        std::getline(input_batch_file, s10,',');
        std::getline(input_batch_file, s11,',');
        std::getline(input_batch_file, s12,',');
        std::getline(input_batch_file, s13);
        //<new>

        batch_dataframe[row].simulationId = s1; 
        batch_dataframe[row].meshFile = s2;
        batch_dataframe[row].boundaryFile = s3;
        batch_dataframe[row].loadFile = s4;
        batch_dataframe[row].paramsFile = s5;
        //batch_dataframe[row].simulationFlag = stoi(s6);
        //batch_dataframe[row].simulationConvergence = stoi(s7);
        //batch_dataframe[row].simulationTime = stold(s8);
        //batch_dataframe[row].simulationFinalTopologyFile = s9;
        //<new>
        batch_dataframe[row].isTopGradSimulation = stoi(s6);
        batch_dataframe[row].inputDensity = s7;
        batch_dataframe[row].simulationFlag = stoi(s8);
        batch_dataframe[row].simulationConvergence = stoi(s9);
        batch_dataframe[row].simulationTime = stold(s10);
        batch_dataframe[row].numberIteration = stoi(s11);
        batch_dataframe[row].finalCompliance = stold(s12);
        //batch_dataframe[row].finalVolumeFraction = stold(s13);
        //batch_dataframe[row].finalLsf = stold(s14);
        batch_dataframe[row].simulationFinalTopologyFile = s13;
        //<new>
    }
    input_batch_file.close();
}

//==================================================================================================
//public member functions
//==================================================================================================
std::string* batch::get_lset_params(int row){
    std::ifstream input_params_file;
    input_params_file.open(params_root + batch_dataframe[row].paramsFile);

    for (int i = 0; i < 15; ++i){
        input_params_file >> lset_params[i];
    }

    input_params_file.close();

    return lset_params;
}

//void batch::update_batch_file(int row, 
//                              int simulationConvergence, 
//                              double simulationTime, 
//                              std::string final_topology_filename){

//<new>
void batch::update_batch_file(int row, 
                              int simulationConvergence, 
                              double simulationTime, 
                              std::string final_topology_filename,
                              int num_iteration,
                              double final_compliance){

    //update batch_dataframe
    batch_dataframe[row].simulationFlag = 1;
    batch_dataframe[row].simulationConvergence = simulationConvergence;
    batch_dataframe[row].simulationTime = simulationTime;
    batch_dataframe[row].simulationFinalTopologyFile = final_topology_filename;
    batch_dataframe[row].finalCompliance = final_compliance;
    batch_dataframe[row].numberIteration = num_iteration;
    
    //update batch csv file
    std::ofstream output_batch_file;
    output_batch_file.open(batch_path);
    for (int i = 0; i < batch_size-1; ++i){
	    
        output_batch_file << batch_dataframe[i].simulationId << ",";
		output_batch_file << batch_dataframe[i].meshFile << ",";
        output_batch_file << batch_dataframe[i].boundaryFile << ",";
        output_batch_file << batch_dataframe[i].loadFile << ",";
        output_batch_file << batch_dataframe[i].paramsFile << ",";
        //output_batch_file << batch_dataframe[i].simulationFlag << ",";
        //output_batch_file << batch_dataframe[i].simulationConvergence << ",";
        //output_batch_file << batch_dataframe[i].simulationTime << ",";
        //output_batch_file << batch_dataframe[i].simulationFinalTopologyFile << "\n";
        //<new>
        output_batch_file << batch_dataframe[i].isTopGradSimulation << ",";
        output_batch_file << batch_dataframe[i].inputDensity << ",";
        output_batch_file << batch_dataframe[i].simulationFlag << ",";
        output_batch_file << batch_dataframe[i].simulationConvergence << ",";
        output_batch_file << batch_dataframe[i].simulationTime << ",";
        output_batch_file << batch_dataframe[i].numberIteration << ",";
        output_batch_file << batch_dataframe[i].finalCompliance << ",";
        output_batch_file << batch_dataframe[i].simulationFinalTopologyFile << "\n";
        //<new>
        
    }
    output_batch_file << batch_dataframe[batch_size-1].simulationId << ",";
	output_batch_file << batch_dataframe[batch_size-1].meshFile << ",";
    output_batch_file << batch_dataframe[batch_size-1].boundaryFile << ",";
    output_batch_file << batch_dataframe[batch_size-1].loadFile << ",";
    output_batch_file << batch_dataframe[batch_size-1].paramsFile << ",";
	//output_batch_file << batch_dataframe[batch_size-1].simulationFlag << ",";
    //output_batch_file << batch_dataframe[batch_size-1].simulationConvergence << ",";
    //output_batch_file << batch_dataframe[batch_size-1].simulationTime << ",";
    //output_batch_file << batch_dataframe[batch_size-1].simulationFinalTopologyFile;
    //<new>
    output_batch_file << batch_dataframe[batch_size-1].isTopGradSimulation << ",";
    output_batch_file << batch_dataframe[batch_size-1].inputDensity << ",";
    output_batch_file << batch_dataframe[batch_size-1].simulationFlag << ",";
    output_batch_file << batch_dataframe[batch_size-1].simulationConvergence << ",";
    output_batch_file << batch_dataframe[batch_size-1].simulationTime << ",";
    output_batch_file << batch_dataframe[batch_size-1].numberIteration << ",";
    output_batch_file << batch_dataframe[batch_size-1].finalCompliance << ",";
    output_batch_file << batch_dataframe[batch_size-1].simulationFinalTopologyFile << "\n";
    //<new>

    output_batch_file.close();

}

void batch::set_simulationFinalTopologyFile(double* density_file, int numElement, std::string file_name){

    std::ofstream final_topology_file;
    final_topology_file.open(final_topology_root + file_name);

    // std::cout << final_topology_root + file_name << std::endl;

    for (int col=0; col<numElement-1; ++col){
        final_topology_file << density_file[col] << ",";
    }
    final_topology_file << density_file[numElement-1] << "\n";
    final_topology_file.close();
}

 //<new>
void batch::set_simulationInitialElementCompliance(double* compliance_file, int numElement, std::string file_name){

    std::ofstream initial_compliance_file;
    initial_compliance_file.open(initial_compliance_root + file_name);

    for (int col=0; col<numElement-1; ++col){
        initial_compliance_file << compliance_file[col] << ",";
    }
    initial_compliance_file << compliance_file[numElement-1] << "\n";
    initial_compliance_file.close();
}

void batch::set_simulationInitialElementTopGrad(double* topgrad_file, int numElement, std::string file_name){

    std::ofstream initial_topgrad_file;
    initial_topgrad_file.open(initial_topgrad_root + file_name);

    for (int col=0; col<numElement-1; ++col){
        initial_topgrad_file << topgrad_file[col] << ",";
    }
    initial_topgrad_file << topgrad_file[numElement-1] << "\n";
    initial_topgrad_file.close();
}
//<new>



std::string batch::get_inputDensity(int row){

    return final_topology_root + batch_dataframe[row].inputDensity;

}


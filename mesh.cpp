#include "mesh.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include <cmath>

mesh::mesh(std::string mesh_file, std::string bdry_file, std::string load_file):
	
	numNode(0), 
	numElement(0), 
	numBoundary(0),
	numLoad(0), 
	nodes(nullptr),
	elements(nullptr), 
	boundaries(nullptr), 
	loads(nullptr),
	nelx(0),
    nely(0),
    elementLengthX(0),
    elementLengthY(0),
	meshFilePath("no file"){
	
	std::ifstream input_mesh_file;
	input_mesh_file.open(mesh_file);	
	std::string stmp;

	if (input_mesh_file.good()){
		input_mesh_file >> numNode;				
		nodes = new nodeStruct [numNode]();
		for (int row = 0; row < numNode; ++row){
	 		input_mesh_file >> nodes[row].nodeId;
	 		input_mesh_file >> nodes[row].xCoord;
	 		input_mesh_file >> nodes[row].yCoord;
		}
		input_mesh_file >> numElement;

		elements = new int [numElement*8]();
		for (int row = 0; row < numElement; ++row){
			for (int col = 0; col < 8; ++col){
				input_mesh_file >> elements[row*8+col]; 
			}
		}
		input_mesh_file.close(); 
	}else{
	 	std::cout << "error: check mesh_file" << std::endl;
	}
	
	// if (input_mesh_file.good()){
	// 	while (!input_mesh_file.eof()){
	// 		getline(input_mesh_file,stmp);
	// 		if (stmp == "$Nodes"){
	// 			input_mesh_file >> numNode;
	// 			nodes = new nodeStruct [numNode]();
	// 			for (int row = 0; row < numNode; ++row){
	// 				input_mesh_file >> nodes[row].nodeId;
	// 				input_mesh_file >> nodes[row].xCoord;
	// 				input_mesh_file >> nodes[row].yCoord;
	// 			}				
			
	// 		} else if (stmp == "$Elements"){
	// 			input_mesh_file >> numElement;
	// 			elements = new int [numElement*8]();
	// 			for (int row = 0; row < numElement; ++row){
	// 				for (int col = 0; col < 8; ++col){
	// 					input_mesh_file >> elements[row*8+col]; 
	// 				}
	// 			}
	// 		}
	// 	}
	// 	input_mesh_file.close();
	// }else{
	// 	std::cout << "error: check mesh_file" << std::endl;
	// }
	
	input_mesh_file.close();

	std::ifstream input_bdry_file;
	input_bdry_file.open(bdry_file);	
	

	if (input_bdry_file.good()){
		input_bdry_file >> numBoundary;				
		boundaries = new boundaryStruct [numBoundary]();
		for (int row = 0; row < numBoundary; ++row){
			input_bdry_file >> boundaries[row].nodeId;
			input_bdry_file >> boundaries[row].nodeDof;
			input_bdry_file >> boundaries[row].boundaryModule;
		}
		input_bdry_file.close(); 
	}else{
		std::cout << "error: check bdry_file" << std::endl;
	}

	std::ifstream input_load_file;
	input_load_file.open(load_file);	
		
	if (input_load_file.good()){
		input_load_file >> numLoad;
		loads = new loadStruct [numLoad]();
		for (int row = 0; row < numLoad; ++row){
			input_load_file >> loads[row].nodeId;
			input_load_file >> loads[row].nodeDof;
			input_load_file >> loads[row].loadModule;
		}
		input_load_file.close(); 
	}else{
		std::cout << "error: check load_file" << std::endl;
	}

	set_numberOfElementsAndLengthPerDiretion();
    meshFilePath = mesh_file;

}

mesh::~mesh(){

	delete [] nodes;
	delete [] elements;	
	delete [] boundaries;
	delete [] loads;	

}
//==================================================================================================
//private member functions
//==================================================================================================

void mesh::set_numberOfElementsAndLengthPerDiretion(){
	elementLengthX = this->nodes[1].xCoord-this->nodes[0].xCoord;
	
	for (int i = 1; i<numNode; ++i){		
		if (nodes[i].xCoord < nodes[i-1].xCoord){
			elementLengthY = this->nodes[i].yCoord - this->nodes[0].yCoord;
			break;
		}
		nelx++;
	}
	nely = numElement/nelx;    
}

//==================================================================================================
// public member functions
//==================================================================================================

void mesh::log_nodes(){
	std::ofstream csv_file_nodes;
	csv_file_nodes.open("result-repo/nodes_mesh.csv");

	csv_file_nodes << "nodeId" << "," << "x" << "," << "y" << "\n";
	for (int row = 0; row < numNode; ++row){
		csv_file_nodes << nodes[row].nodeId << ",";
		csv_file_nodes << nodes[row].xCoord << ",";
		csv_file_nodes << nodes[row].yCoord << "\n";
	}
	csv_file_nodes.close();
}

void mesh::log_elements(){
	std::ofstream csv_file_elements;
	csv_file_elements.open("result-repo/elements_mesh.csv");

	for (int row = 0; row < numElement; ++row){
		for (int col = 0; col < 7; ++col){
			csv_file_elements << elements[row*8+col] << ",";
		}
		csv_file_elements << elements[row*8+7] << "\n";
	}
	csv_file_elements.close();
}

void mesh::log_boundaries(){
	std::ofstream csv_file_boundaries;
	csv_file_boundaries.open("result-repo/boundaries_mesh.csv");

	csv_file_boundaries << "nodeId" << "," << "nodeDof" << "," << "boundaryModule" << "\n";
	for (int row = 0; row < numBoundary; ++row){
		csv_file_boundaries << boundaries[row].nodeId << ",";
		csv_file_boundaries << boundaries[row].nodeDof << ",";
		csv_file_boundaries << boundaries[row].boundaryModule << "\n";
	}
	csv_file_boundaries.close();
}

void mesh::log_loads(){
	std::ofstream csv_file_loads;
	csv_file_loads.open("result-repo/loads_mesh.csv");

	csv_file_loads << "nodeId" << "," << "nodeDof" << "," << "loadModule" << "\n";
	for (int row = 0; row < numLoad; ++row){
		csv_file_loads << loads[row].nodeId << ",";
		csv_file_loads << loads[row].nodeDof << ",";
		csv_file_loads << loads[row].loadModule << "\n";
	}
	csv_file_loads.close();
}
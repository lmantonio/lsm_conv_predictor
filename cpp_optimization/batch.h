#ifndef BATCH_H
#define BATCH_H

#include <string>

struct batchData{
	std::string simulationId;
    std::string meshFile;
	std::string boundaryFile;
    std::string loadFile;
    std::string paramsFile;
    // <new> 
    int isTopGradSimulation;
    std::string inputDensity;
    // </new>
    int simulationFlag;
    int simulationConvergence;
    double simulationTime;
    // <new> 
    int numberIteration;
    double finalCompliance;
    // </new>
    std::string simulationFinalTopologyFile;

};

class batch{
    private:
        
        int batch_size;
        std::string batch_root;
        std::string batch_path;
        std::string mesh_root;
        std::string bdry_root;
        std::string load_root;
        std::string params_root;
        std::string final_topology_root;
        std::string initial_compliance_root;
        std::string initial_topgrad_root;
        std::string *lset_params;
        batchData *batch_dataframe;

        void set_batch_size();
        void set_batch_dataframe();
    
    public:
        batch(std::string batch_file_path, std::string batch_name);
        ~batch();

        inline int get_batch_size() const{return batch_size;};
        inline std::string get_simulationId(int row) {return batch_dataframe[row].simulationId;};
        inline std::string get_meshFilePath(int row) {return (mesh_root + batch_dataframe[row].meshFile);};
        inline std::string get_boundaryFilePath(int row) {return (bdry_root + batch_dataframe[row].boundaryFile);};
        inline std::string get_loadFilePath(int row) {return (load_root + batch_dataframe[row].loadFile);};
        inline std::string get_paramsFilePath(int row) {return (params_root + batch_dataframe[row].paramsFile);}; 
        inline int get_simulationFlag(int row) {return batch_dataframe[row].simulationFlag;};
       
        std::string* get_lset_params(int row);
        void set_simulationFinalTopologyFile(double* density_file, int numElement, std::string file_name);
        void set_simulationInitialElementCompliance(double* compliance_file, int numElement, std::string file_name);
        void set_simulationInitialElementTopGrad(double* topgrad_file, int numElement, std::string file_name);
        //void update_batch_file(int row, int simulationConvergence, double simulationTime, std::string final_topology_filename);
        //<new>
        void update_batch_file(int row, 
                               int simulationConvergence, 
                               double simulationTime, 
                               std::string final_topology_filename,
                               int num_iteration,
                               double final_compliance);
        inline int get_lsetMode(int row) {return batch_dataframe[row].isTopGradSimulation;};
        std::string get_inputDensity(int row);
        //</new> 
};

#endif
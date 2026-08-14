#ifndef MESH_H
#define MESH_H

#include <string>

struct nodeStruct{
	int nodeId;
	double xCoord;
	double yCoord;
};

struct boundaryStruct{
	int nodeId;
	int nodeDof;
	double boundaryModule;
};

struct loadStruct{
	int nodeId;
	int nodeDof;
	double loadModule; 
};

class mesh{
	private:
		int numNode;
        int numElement;
        int numBoundary;
        int numLoad;
		nodeStruct *nodes;
		int *elements;
		boundaryStruct *boundaries;
		loadStruct *loads;       
        int nelx;
        int nely;
        double elementLengthX;
        double elementLengthY; 
        std::string meshFilePath;

        void set_numberOfElementsAndLengthPerDiretion();
       

	public:
    	mesh(std::string mesh_file, std::string bdry_file, std::string load_file);
    	~mesh();
        
        inline const int get_numNode() const{return numNode;};
        inline const int get_numElement() const{return numElement;};
        inline const int get_numBoundary() const{return numBoundary;};
        inline const int get_numLoad() const{return numLoad;};
    	inline const nodeStruct* get_Nodes() const{return nodes;};
    	inline const int* get_Elements() const{return elements;};
    	inline const boundaryStruct* get_Boundaries() const{return boundaries;};
    	inline const loadStruct* get_Loads() const{return loads;};
        inline const int get_nelx() const{return nelx;};
        inline const int get_nely() const{return nely;};
        inline const double get_elementLengthX() const{return elementLengthX;};
        inline const double get_elementLengthY() const{return elementLengthY;};
        inline const std::string get_meshFilePath() const{return meshFilePath;};

        void log_nodes();
        void log_elements();
        void log_boundaries();
        void log_loads();
                

};
#endif
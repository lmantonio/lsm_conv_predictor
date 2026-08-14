# LSM Convergence Predictor

This repository contains the source code and implementation details associated with the PhD research:

> **A Proposal for Reducing Dependency on Parameter Selection in Structural Topology Optimization by Level Set Method**  
> **Author:** Leonardo Machado Antonio

## Overview
This repository contains a collection of codes divided into two main modules:

1. **Level Set Topology Optimization:** C++ implementation (header and source files) designed to execute 2D level set topology optimization on structures.
2. **Machine Learning Convergence Prediction Training:** Python implementation for training machine learning models to predict the convergence or non-convergence of 2D level set topology optimizations.

---

## Requirements
* C++ compiler (supporting standard C++ headers)
* Python environment with standard data science and machine learning libraries (check specific imports in the scripts)

---

## Usage

### 1. Level Set Topology Optimization (C++)
1. Compile the source files along with a main file (`topop.cpp`) to generate the executable.
2. In the root directory, create the following folder structure: `repo`, `meshes`, `loads`, `params`, `boundaries`, `initial_compliance`, `initial_toplgrad`, and `final_topologies`. Place the compiled executable in the root directory as well.
3. Place a CSV file containing the simulation batch (without headers) inside the `repo` folder.
4. Run the executable by passing the CSV file name as an argument to start batch processing:
   ```bash
   ./your_executable batch_file.csv


### 2. Machine Learning Training (Python)
1. Place the `gen2_repo` folder (containing the processed dataset) in the same root directory as the Python scripts.
2. Run the Python script in your terminal or command prompt, including the arguments for the dataset path (`gen2_repo`) and the dataset fragment ID.

**Dataset Naming Convention:**
* Files prefixed with `ml` correspond to classical machine learning algorithms.
* Files prefixed with `dl` correspond to Deep Neural Networks (DNN) and Transformer algorithms.

Example format: `x_train_gen2_f_c_dnn_global.csv`
* **f/p:** Full feature set (`f`) or partial feature set (`p`).
* **c/t:** Compliance field only (`c`) or compliance and topological gradient data (`t`).
* **global:** Balancing/stratification type.
* **x/y:** Input parameters/training dataset (`x`) or corresponding convergence flag (`y`).

**Outputs:**
Upon completion, the training script saves its results in a `results` folder containing the model name and the following subfolders/files:
* **Best Results (`best_result/`):** 
  * `best_model.pkl` (trained model object)
  * Best parameters CSV
  * Threshold vs. performance metric data
* **General Results (`general_results/`):** 
  * Precision-recall curve data
  * ROC curve data
  * Threshold evaluation results data



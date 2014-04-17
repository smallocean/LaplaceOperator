#pragma once
#include "DS.h"
#include "PCL.h"


//生成2维的稀疏矩阵
void generate_laplace_matrix_sparse_matlab_dim2(PCloud& pcloud, double h, double rho, vector<unsigned int>& II, vector<unsigned int>& JJ, vector<double>& SS);

//生成3维的稀疏矩阵
void generate_laplace_matrix_sparse_matlab_dim3(PCloud& pcloud, double h, double rho, vector<unsigned int>& II, vector<unsigned int>& JJ, vector<double>& SS);

//生成k维的稀疏矩阵
void generate_laplace_matrix_sparse_matlab_dimk(PCloud& pcloud, double h, double rho, unsigned int tdim, vector<unsigned int>& II, vector<unsigned int>& JJ, vector<double>& SS);

////生成图Laplace矩阵
void generate_graph_laplace_matrix_sparse_matlab(PCloud& pcloud, double h, double rho, unsigned int tdim, vector<unsigned int>& II, vector<unsigned int>& JJ, vector<double>& SS);

//生成测地线距离的Laplace矩阵
void generate_arbdist_graph_laplace_matrix_sparse_matlab(PCloud& pcloud, double h, double rho, unsigned int tdim, vector<unsigned int>& II, vector<unsigned int>& JJ, vector<double>& SS);

//生成核函数矩阵
void generate_kernel_matrix_sparse_matlab(PCloud& pcloud, double h, double rho, unsigned int tdim, vector<unsigned int>& II, vector<unsigned int>& JJ, vector<double>& SS);


//估计切平面
void estimateTangentSpace(const dPoint& pt, const vector<dPoint>& neighbor_pts, double h, unsigned int tdim, vector<dVector>& tspace);

//单一标记量
double simplex_sign_volume(vector<dPoint> points);



//生成点云Laplace矩阵
void generate_pcdlaplace_symmetric_matrix_sparse(
        PCloud& pcloud,
        double h,
        double rho,
        vector<unsigned int>& II,
        vector<unsigned int>& JJ,
        vector<double>& SS,
        vector<double>& BB);

//计算矩阵QB
void compute_QB(PCloud& pcloud, double h, double rho, vector<unsigned int>& II, vector<unsigned int>& JJ,vector<double>& SS,vector<double>& BB);



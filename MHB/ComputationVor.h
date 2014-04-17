#pragma once
#include "DS.h"
#include "PCL.h"


//����2ά��ϡ�����
void generate_laplace_matrix_sparse_matlab_dim2(PCloud& pcloud, double h, double rho, vector<unsigned int>& II, vector<unsigned int>& JJ, vector<double>& SS);

//����3ά��ϡ�����
void generate_laplace_matrix_sparse_matlab_dim3(PCloud& pcloud, double h, double rho, vector<unsigned int>& II, vector<unsigned int>& JJ, vector<double>& SS);

//����kά��ϡ�����
void generate_laplace_matrix_sparse_matlab_dimk(PCloud& pcloud, double h, double rho, unsigned int tdim, vector<unsigned int>& II, vector<unsigned int>& JJ, vector<double>& SS);

////����ͼLaplace����
void generate_graph_laplace_matrix_sparse_matlab(PCloud& pcloud, double h, double rho, unsigned int tdim, vector<unsigned int>& II, vector<unsigned int>& JJ, vector<double>& SS);

//���ɲ���߾����Laplace����
void generate_arbdist_graph_laplace_matrix_sparse_matlab(PCloud& pcloud, double h, double rho, unsigned int tdim, vector<unsigned int>& II, vector<unsigned int>& JJ, vector<double>& SS);

//���ɺ˺�������
void generate_kernel_matrix_sparse_matlab(PCloud& pcloud, double h, double rho, unsigned int tdim, vector<unsigned int>& II, vector<unsigned int>& JJ, vector<double>& SS);


//������ƽ��
void estimateTangentSpace(const dPoint& pt, const vector<dPoint>& neighbor_pts, double h, unsigned int tdim, vector<dVector>& tspace);

//��һ�����
double simplex_sign_volume(vector<dPoint> points);



//���ɵ���Laplace����
void generate_pcdlaplace_symmetric_matrix_sparse(
        PCloud& pcloud,
        double h,
        double rho,
        vector<unsigned int>& II,
        vector<unsigned int>& JJ,
        vector<double>& SS,
        vector<double>& BB);

//�������QB
void compute_QB(PCloud& pcloud, double h, double rho, vector<unsigned int>& II, vector<unsigned int>& JJ,vector<double>& SS,vector<double>& BB);



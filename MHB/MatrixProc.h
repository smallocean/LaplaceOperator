#pragma once
#include <vector>
class MatrixProc
{
public:
	MatrixProc(void);
	~MatrixProc(void);
};

int convertIJS2colMatrix(std::vector<unsigned int>, std::vector<unsigned int> JJ, std::vector<double> SS, double* a, int m, int n);

int convertVector2Array(std::vector<double> v, double* y);

std::vector<double> readDoubleArray(const char* filename);

bool writeDoubleArray(std::vector<double>& param, const char* filename);

//����Voronoiͼ���
double voronoi_cell_area(std::vector<std::vector<double> > &m_points_2d);

int readPurePointsFile(const char* filename, std::vector<std::vector<double> > &m_points);

void cal_pcd_laplace(std::vector<std::vector<double> > &m_points,
        std::vector<unsigned int> &m_II,
        std::vector<unsigned int> &m_JJ,
        std::vector<double> &m_SS);



void cal_sym_pcd_laplace(std::vector<std::vector<double> > &m_points,
        std::vector<unsigned int> &m_II,
        std::vector<unsigned int> &m_JJ,
        std::vector<double> &m_SS,
        std::vector<double> &m_BB);


bool save_IJS(
        const char* filename,  //�ļ�·��
        std::vector<unsigned int>& II, //�����к�
        std::vector<unsigned int>& JJ, //�����к�
        std::vector<double>& SS);     //����Ȩֵ ��1/(4pi*t^2)*exp(-d/(4t^2))

int solve_eigen(double* v,
        double* lambda, long int n);

int convert_IJS_2_col_sq_matrix(
        double* result,
        std::vector<unsigned int> II,
        std::vector<unsigned int> JJ,
        std::vector<double> SS, int n );

int writeDoubleRowSqMatrix(double* p, int n, const char* filename);

int writeDoubleRow(double* p, int n, const char* filename);

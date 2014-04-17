

#include <stdlib.h>

#include <stdio.h>


#include <string>
#include <fstream>
#include <iostream>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>
#include "MatrixProc.h"
#include "ComputationVor.h"
#include "DS.h"
#include "PCL.h"
#include "IOProc.h"

using namespace std;

/*
 *
 */


void cal_sym_pcd_qb(std::vector<std::vector<double> > &m_points, std::vector<unsigned int> &m_II, std::vector<unsigned int> &m_JJ, std::vector<double> &m_SS, std::vector<double> &m_BB) {


	//m_II
	//m_JJ
	//m_SS
	//M_BB 保存边界点的面积
    if (m_points.size() > 10) {
        double x, y, z;

        PCloud *pc;
        int i, size;
        double* ar;
        double avgs;

        m_II.resize(0);
        m_JJ.resize(0);
        m_SS.resize(0);
        m_BB.resize(0);

        size = m_points.size();  //点的个数
        ar = new double[3 * size];//点坐标

        i = 0;

        std::vector<double> pt;

        BOOST_FOREACH(pt, m_points) {
            x = pt[0];
            y = pt[1];
            z = pt[2];

            ar[i] = x;
            ar[i + size] = y;
            ar[i + size * 2] = z;

            i++;
        }
        pc = new PCloud(ar, size, 3);
        avgs = pc->average_size(10);

        cerr << "cal_sym_pcd_qb to reserve for IJSB .\n";
        m_II.reserve(m_points.size()*200);
        cerr << "reservation for I done " << m_points.size()*180 << endl;
        m_JJ.reserve(m_points.size()*200);
        cerr << "reservation for J done " << m_points.size()*180 << endl;
        m_SS.reserve(m_points.size()*200);
        cerr << "reservation for S done " << m_points.size()*180 << endl;
        m_BB.reserve(m_points.size());
        cerr << "reservation for B done " << m_points.size() << endl;
        cerr << "cal_sym_pcd_qb reservation done.\n";

        compute_QB(*pc, avgs * 2, 3,m_II, m_JJ, m_SS, m_BB);

        /* put average area for all boundry vertices */
        double val, sum = 0;
        i = 0;

        BOOST_FOREACH(val, m_BB) {
            if (val > 0) {
                sum += val;
                ++i;
            } else
                printf("cal_sym_pcd_qb: boundry detected.\n");
        }
        sum /= i;

        for (i = 0; i < m_BB.size(); ++i) {
            if (m_BB[i] <= 0) m_BB[i] = sum;
        }

		
        /* add area factor for all non-diag elements
         * currently assume that
         * PCDLap code generate symmetric connectivity
         */
        std::vector<double> diag;
        diag.resize(m_points.size(), 0);
        for (i = 0; i < m_II.size(); ++i) {
            if (m_II[i] == m_JJ[i]) {
                cerr << "cal_sym_pcd_qb: errornous diagonal element \n";
                continue;
            }

            if (m_JJ[i] > m_II[i]) {
                m_SS[i] *= m_BB[m_JJ[i] - 1];
                m_SS[i] *= m_BB[m_II[i] - 1];
            } else {
                m_SS[i] *= m_BB[m_II[i] - 1];
                m_SS[i] *= m_BB[m_JJ[i] - 1];
            }
            diag[m_II[i] - 1] += m_SS[i];
        }
        /* add the diagonal elements */
        for (i = 0; i < diag.size(); ++i) {
            m_II.push_back(i + 1);
            m_JJ.push_back(i + 1);
            m_SS.push_back(0 - diag[i]);
        }

        delete pc;
        delete[] ar;
    }
}


int main(int argc, char** argv) {

    //if (argc < 4) {
    //    printf("not enough parameters: %d\n", argc);
    //    printf("Usage: pcdqb.exe MModel(input) QFile(output) BFile(output)");
    //    return 0;
    //}

    pcdwrapper pw;
	//cout<<"please input the model(.m):"<<endl;
	//string input,outputQ,outputB;
	//cin>>input;
	//char * p =input;
   // pw.readMFile(argv[1]);
	string input = "D:\Laurana50k.m";
	char *a = "D:\source\\f.m";
	char *b = "D:\\source\\Q1.txt";
	char *c = "D:\\source\\B1.txt";
	pw.readMFile(a);
    printf("%d vertices read.\n", pw.m_points.size());
    if (pw.m_points.size() <= 0) {
        printf("Reading Error\n");
    }

    std::vector<unsigned int> I, J;
    std::vector<double> S, B;

    I.clear();
    J.clear();
    S.clear();
    B.clear();

    cal_sym_pcd_qb(pw.m_points, I, J, S, B);
	//cout<<"please input the route of outputB"<<endl;
	//cin>>outputB;
    //    save_IJS("pcdlaplace_vd_matrix.txt", I, J, S) ? printf("IJSB true\n") : printf("IJSB false\n");
 //   printf("saving q matrix to: %s\n", argv[2]);
	printf("saving q matrix to: %s\n", b);
    save_IJS(argv[2], I, J, S) ? printf("save sym IJS true\n") : printf("save sym IJS false\n");
    //    writeDoubleArray(B, "pcdlaplace_vd_mass.txt");
   // printf("saving b matrix to: %s\n", argv[3]);
	printf("saving b matrix to: %s\n", c);
    writeDoubleArray(B, argv[3]);

    int np;
    np = pw.m_points.size();

    double *p = NULL, *d = NULL;
    /*
    p = new double[np*np];
    d = new double[np];

    convert_IJS_2_col_sq_matrix(p, I, J, S, np);
    solve_eigen(p, d, np);

    printf("saving eigen vectors to: %s\n", argv[4]);
    writeDoubleRowSqMatrix(p, np, argv[4]);
    printf("saving eigen values to: %s\n", argv[5]);
    writeDoubleRow(d, np, argv[5]);

    delete[] p;
    delete[] d;
     */
    return (EXIT_SUCCESS);
}


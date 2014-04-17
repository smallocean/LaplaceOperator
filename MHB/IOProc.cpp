#include "IOProc.h"


#include "MatrixProc.h"
#include "ComputationVor.h"

#include <vector>
#include <fstream>
#include <iostream>
#include <string>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>


extern "C" { int dgemv_(char *, long int *, long int *,double *, double *, long int *, double *, long int *, double *, double *, long int *);
}

pcdwrapper::pcdwrapper(void) {
}

pcdwrapper::~pcdwrapper(void) {
}

bool pcdwrapper::readMFile(const char *filename) {
    std::ifstream ifs;
    std::string line;

    std::string s;
    std::vector<std::string> stvec;

    std::vector<std::string> filtervec;

    std::vector<double> pt;

    ifs.open(filename);

    double x, y, z, val;


    if (ifs.is_open()) {
        m_points.clear();
        pt.clear();

        do {
            std::getline(ifs, line);
            if (line.size() <= 0) continue;
            boost::algorithm::split(stvec, line, boost::algorithm::is_any_of(" \n"));
//
//            if (stvec.size() < 5) {
//                continue;
//            }

            filtervec.clear();
            BOOST_FOREACH(s, stvec) {
                if(s!="") {
                    filtervec.push_back(s);
                }
            }

            if(filtervec.size()<5) {
                continue;
            }

            try {
                if (filtervec[0] == "Vertex") {
                    x = boost::lexical_cast<double>(filtervec[2]);
                    y = boost::lexical_cast<double>(filtervec[3]);
                    z = boost::lexical_cast<double>(filtervec[4]);

                    /*
                    if (stvec.size() > 5) {
                        val = boost::lexical_cast<double>(stvec[5]);
                    }
                     */

                    pt.clear();

                    pt.push_back(x);
                    pt.push_back(y);
                    pt.push_back(z);
                    //pt.push_back(val);

                    m_points.push_back(pt);
					
	            }
            } catch (boost::bad_lexical_cast&) {
                printf("pcdwrapper::readMFile: lexical_cast fail: %d: %s %s %s\n",
                        stvec.size(),
                        stvec[2].c_str(), stvec[3].c_str(), stvec[4].c_str());
            }
        } while (!ifs.eof());
    }
    
	if (m_points.size() <= 0) {
        return false;
    }

    return true;
}

bool pcdwrapper::readMatFile(const char* filename) {
    std::ifstream ifs;
    std::string line;

    std::string s;
    std::vector<std::string> stvec;

    ifs.open(filename);

    double val;


    if (ifs.is_open()) {
        m_val.clear();

        std::getline(ifs, line);
        if (line.size() <= 0) return false;

        boost::algorithm::split(stvec, line, boost::algorithm::is_any_of(" \n"));

        BOOST_FOREACH(s, stvec) {
            try {
                val = boost::lexical_cast<double>(s);
                m_val.push_back(val);
            } catch (boost::bad_lexical_cast&) {
                printf("pcdwrapper::readMatFile: lexical_cast fail at No. %d\n", m_val.size());
            }
        }
    } else {
        return false;
    }
    if (m_val.size() <= 0) {
        return false;
    }

    return true;
}

std::vector<double> pcdwrapper::multiplyMatrix() {
    std::vector<double> result;

    if (m_II.size() <= 0) {
        return result;
    }

    long int size = m_points.size(), inc = 1;
    double *A = new double[size * size];
    double *xarray = new double[size];
    double *yarray = new double[size];
    double alpha = 1, beta = 0;
    int i, j;
    std::vector<double> rt;
    int lapack_result;

    memset(A, 0, sizeof (double) * size * size);

    convertIJS2colMatrix(m_II, m_JJ, m_SS, A, size, size);
    convertVector2Array(m_val, xarray);

    lapack_result = dgemv_("N",
            &size, //m
            &size, //n
            &alpha, // alpha
            A,
            &size, //lda
            xarray, // x
            &inc, // incx
            &beta, //beta
            yarray, // y
            &inc); //incy

    for (i = 0; i < size; i++) {
        result.push_back(yarray[i]);
    }

    delete[] A;
    delete[] xarray;
    delete[] yarray;

    return result;
}

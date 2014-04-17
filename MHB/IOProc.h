#pragma once

#include <vector>

class pcdwrapper {
public:
    pcdwrapper(void);
    ~pcdwrapper(void);

public:
	//读取.m文件 m文件格式  vertex           face
    bool readMFile(const char *filename);
    bool readMatFile(const char *filename);

	//计算点云
    void calculate_pcd();  

	//保存点云数据
    std::vector<std::vector<double> > m_points;
    
	//
	std::vector<unsigned int> m_II, m_JJ;
    std::vector<double> m_SS;
    std::vector<double> m_val;

    std::vector<std::vector<unsigned int> > m_sorted_piv;

    std::vector<double> multiplyMatrix();
};

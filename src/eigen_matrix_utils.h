#ifndef EIGEN_MATRIX_UTILS_INCLUDED
#define EIGEN_MATRIX_UTILS_INCLUDED

#include <eigen3/Eigen/Dense>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/split.hpp>
#include <exception>

#define MAXBUFSIZE  ((int) 1e6)




typedef std::vector<std::string> string_vector;


void removeRow(Eigen::MatrixXd& matrix, unsigned int rowToRemove)
{
	unsigned int numRows = matrix.rows()-1;
	unsigned int numCols = matrix.cols();
	
	if( rowToRemove < numRows )
		matrix.block(rowToRemove,0,numRows-rowToRemove,numCols) = matrix.block(rowToRemove+1,0,numRows-rowToRemove,numCols);
	
	matrix.conservativeResize(numRows,numCols);
}

void removeColumn(Eigen::MatrixXd& matrix, unsigned int colToRemove)
{
	unsigned int numRows = matrix.rows();
	unsigned int numCols = matrix.cols()-1;
	
	if( colToRemove < numCols )
		matrix.block(0,colToRemove,numRows,numCols-colToRemove) = matrix.block(0,colToRemove+1,numRows,numCols-colToRemove);
	
	matrix.conservativeResize(numRows,numCols);
}



Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> readDynamicMatrix(std::string filename)
{
    int cols = 0, rows = 0;
    double buff[MAXBUFSIZE];

    // Read numbers from file into buffer.
    std::ifstream infile;
    infile.open(filename);

    if (!infile.is_open()){
        std::cerr << "ERROR: readDynamicMatrix: couldn't open file: " << filename << std::endl;
        throw std::exception();
    }

    while (! infile.eof())
    {
        std::string line;
        getline(infile, line);

        int temp_cols = 0;
        std::stringstream stream(line);
        while(! stream.eof())
            stream >> buff[cols*rows+temp_cols++];

        if (temp_cols == 0)
            continue;

        if (cols == 0)
            cols = temp_cols;

        rows++;
    }

    infile.close();

    rows--;

    // Populate matrix with numbers.
    Eigen::MatrixXd result(rows,cols);
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            result(i,j) = buff[ cols*i+j ];

    return result;
}


template <typename SCALAR,unsigned int DIM_ROW, unsigned int DIM_COL>
Eigen::Matrix<double, DIM_ROW, DIM_COL> readMatrix(std::string filename)
{

    std::ifstream input(filename.c_str(), std::ios::in);

    //const unsigned int dim = 11;

    Eigen::Matrix<SCALAR, DIM_ROW, DIM_COL> mat = Eigen::Matrix<SCALAR, DIM_ROW, DIM_COL>::Zero();
    //mean = Eigen::Matrix<double, dim, 1>::Zero();


    string_vector SplitVec;
    std::string lineStr;
    unsigned long length, lineNum = 0;

    while ( !input.eof() )
    {
        getline(input, lineStr);
        std::string resStr;
        //lineNum++;

        length = lineStr.length();


        if (length > 0)
        {
            boost::split( SplitVec, lineStr, boost::is_any_of("\t") );

            for (unsigned int ii = 0; ii < SplitVec.size(); ii++)
            {
                //cout << lineNum << "\t" << ii << "\t" << SplitVec[ii] << endl;
                mat(lineNum, ii) =   boost::lexical_cast<SCALAR>(SplitVec[ii]);
            }
            lineNum++;
        }
    }

    return mat;
}



#endif // EIGEN_MATRIX_UTILS_INCLUDED

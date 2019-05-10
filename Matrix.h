#pragma once

#include "Vector.h"
#include <iostream>
#include <vector>

template <class T>
class Matrix
{
	private:
		int rows; //num of rows
		int cols; //num of columns
		std::vector<vectort> major_row;

	public:

		Matrix<T>() : rows(1), cols(1)
		{
			major_row = std::vector<vectort>(1, vectort(1, 0));
			//major_row.push_back(vectort(1, 0));
		}

		Matrix<T>(int r, int c) : rows(r), cols(c)
		{
			major_row = std::vector<vectort>(rows, vectort(cols, 0));
		}

		Matrix<T>(const Matrix& m) : rows(m.rows), cols(m.cols), major_row(m.major_row)
		{}

		int get_rows() const
		{
			return rows;
		}

		int get_cols() const
		{
			return cols;
		}

		const T& operator () (int i, int j) const
		{
			return major_row[i][j];
		}

		T& operator () (int i, int j)
		{
			return major_row[i][j];
		}
		
		T at(int i, int j) const
		{
			if (i < rows && j < cols && i >= 0 && j >= 0)
			{
				return major_row[i][j];
			}
			else
			{
				return NULL;
			}
		}

		Matrix<T> transpose() const
		{
			Matrix<T> transed(cols, rows);
			for (int i = 0; i < rows; ++i)
			{
				for(int j = 0; j < cols; ++j)
				{
					transed(j, i) = major_row[i][j];
				}
			}
			return transed;
		}

		void transpose_in_place(Matrix m){}//sudo TODO

		Matrix<T> operator + (const Matrix &m) const
		{
			Matrix<T> sum(rows, cols);
			for (int i = 0; i < rows; ++i)
			{
				for(int j = 0; j < cols; ++j)
				{
					sum(i, j) = major_row[i][j] + m(i, j);
				}
			}
			return sum;	
		}

		Matrix<T> operator * (const Matrix &m) const
		{
			Matrix<T> prod(rows, cols);
			T tmp_val;
			for (int i = 0; i < rows; ++i)
			{
				for(int j = 0; j < m.get_cols(); ++j)
				{
					tmp_val = 0;
					for (int k = 0; k < cols; ++k)
					{
						tmp_val += major_row[i][k] * m(k, j);
					}
					prod(i, j)=tmp_val;
				}
			}
			return prod;		
		}

		Matrix<T> operator - (const Matrix &m) const
		{
			Matrix<T> diff(rows, cols);
			for (int i=0; i<rows; i++)
			{
				for(int j = 0; j < cols; ++j)
				{
					diff(i, j)=major_row[i][j]-m(i, j);
				}
			}
			return diff;
		}

		Matrix<T> operator - () const
		{
			Matrix<T> neg(rows, cols);
			for (int i = 0; i < rows; ++i)
			{
				for(int j = 0; j < cols; ++j)
				{
					neg(i, j) = -major_row[i][j];
				}
			}
			return neg;
		}

		friend Vector<T> operator * (const Matrix<T> &m, const Vector<T> &v)
		{
			Vector<T> res(m.get_rows());
			T tmp_val;
			for (int i = 0; i < m.get_rows(); ++i)
			{
				tmp_val = 0;
				for(int j = 0; j < m.get_cols(); ++j)
				{
					tmp_val += m(i, j) * v(j);
				}
				res(i) = tmp_val;
			}
			return res;	
		}

		friend Matrix<T> operator * (const Matrix& m, const T& s)
		{
			Matrix<T> prod(m.get_rows(), m.get_cols());
			for (int i = 0; i < m.get_rows(); ++i)
			{
				for(int j = 0; j < m.get_cols(); ++j)
				{
					prod(i, j) = m(i, j) * s;
				}
			}
			return prod;
		}

		friend Matrix<T> operator * (const T& s, const Matrix& m)
		{
			return m * s;
		}

};

template <class T>
std::ostream& operator<<(std::ostream& strm, const Matrix<T>& m)
{
	strm << "[";
	for (int i = 0; i < m.get_rows(); ++i)
	{
		strm << "(";
		for (int j = 0; j < m.get_cols() - 1; ++j)
		{
			strm << m(i,j) << ", ";
		}
		strm << m(i, m.get_cols() - 1) << ")\n";
	}
	strm << "]" << std::endl;
	return strm;
}
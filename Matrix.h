#pragma once

#include "Vector.h"
#include <algorithm>
#include <iostream>
#include <vector>

template <class T>
class Matrix {
    using vector_t = std::vector<T>;

private:
    int rows;                        // num of rows
    int cols;                        // num of columns
    std::vector<vector_t> major_row; // the matrix's entries, stored as a vector of rows of the matrix, in row major fashion.

    // max(s: 2^s<=x) x must be nonegative
    int round_up_power_2(int x)
    {
        int i = 0;
        while (x > (1 << i)) {
            ++i;
        }
        return i;
    }

    Matrix<T> mul(const Matrix<T>& A, const Matrix<T>& B, int power)
    {
        int size = 1 << power;
        Matrix<T> prod(size, size);
        if (size == 1) {
            prod(0, 0) = A(0, 0) * B(0, 0);
            return prod;
        }
        int half_size = size / 2;
        Matrix<T> A11(half_size, half_size);
        Matrix<T> A12(half_size, half_size);
        Matrix<T> A21(half_size, half_size);
        Matrix<T> A22(half_size, half_size);
        Matrix<T> B11(half_size, half_size);
        Matrix<T> B12(half_size, half_size);
        Matrix<T> B21(half_size, half_size);
        Matrix<T> B22(half_size, half_size);
        Matrix<T> C11(half_size, half_size);
        Matrix<T> C12(half_size, half_size);
        Matrix<T> C21(half_size, half_size);
        Matrix<T> C22(half_size, half_size);
        Matrix<T> M1(half_size, half_size);
        Matrix<T> M2(half_size, half_size);
        Matrix<T> M3(half_size, half_size);
        Matrix<T> M4(half_size, half_size);
        Matrix<T> M5(half_size, half_size);
        Matrix<T> M6(half_size, half_size);
        Matrix<T> M7(half_size, half_size);

        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                if (i < half_size && j < half_size) {
                    A11(i, j) = A(i, j);
                    B11(i, j) = B(i, j);
                }
                if (i < half_size && j >= half_size) {
                    A12(i, j - half_size) = A(i, j);
                    B12(i, j - half_size) = B(i, j);
                }
                if (i >= half_size && j < half_size) {
                    A21(i - half_size, j) = A(i, j);
                    B21(i - half_size, j) = B(i, j);
                }
                if (i >= half_size && j >= half_size) {
                    A22(i - half_size, j - half_size) = A(i, j);
                    B22(i - half_size, j - half_size) = B(i, j);
                }
            }
        }

        M1 = mul(A11 + A22, B11 + B22, power - 1);
        M2 = mul(A21 + A22, B11, power - 1);
        M3 = mul(A11, B12 - B22, power - 1);
        M4 = mul(A22, B21 - B11, power - 1);
        M5 = mul(A11 + A12, B22, power - 1);
        M6 = mul(A21 - A11, B11 + B12, power - 1);
        M7 = mul(A12 - A22, B21 + B22, power - 1);

        C11 = M1 + M4 - M5 + M7;
        C12 = M3 + M5;
        C21 = M2 + M4;
        C22 = M1 - M2 + M3 + M6;

        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                if (i < half_size && j < half_size) {
                    prod(i, j) = C11(i, j);
                }
                if (i < half_size && j >= half_size) {
                    prod(i, j) = C12(i, j - half_size);
                }
                if (i >= half_size && j < half_size) {
                    prod(i, j) = C21(i - half_size, j);
                }
                if (i >= half_size && j >= half_size) {
                    prod(i, j) = C22(i - half_size, j - half_size);
                }
            }
        }

        return prod;
    }

public:
    // default constructor. Constructs a 1X1 with 0
    explicit Matrix() : rows(1), cols(1), major_row(std::vector<vector_t>(1, vector_t(1, 0)))
    {
    }

    // constructs an rXc matrix, filled with 0's
    explicit Matrix(int r, int c) : rows(r), cols(c), major_row(std::vector<vector_t>(rows, vector_t(cols, 0)))
    {
    }

    // constructor by const reference to std::vector of vector_t (the entries)
    explicit Matrix(const std::vector<vector_t>& entries) : rows(entries.size()), cols(entries[0].size()), major_row(entries)
    {
    }

    // constructor by rvalue reference to std::vector of vector_t (the entries)
    explicit Matrix(std::vector<vector_t>&& entries) : rows(entries.size()), cols(entries[0].size()), major_row(std::move(entries))
    {
    }

    // copy constructor. Constructs a copy of a given matrix
    Matrix(const Matrix& m) : rows(m.rows), cols(m.cols), major_row(m.major_row)
    {
    }

    // move constructor by rvalue reference of a different matrix
    Matrix(Matrix&& m) noexcept : rows(m.rows), cols(m.cols), major_row(std::move(m.major_row))
    {
    }

    // copy assignment by const reference
    Matrix& operator=(const Matrix& m)
    {
        rows = m.rows;
        cols = m.cols;
        major_row = m.major_row;
    }

    // move assignment by rvalue reference
    Matrix& operator=(Matrix&& m) noexcept
    {
        major_row = std::move(m.major_row);
    }

    // return whether or not the given matrix has the same entries.
    bool operator==(const Matrix& m) const
    {
        return major_row == m.major_row;
    }

    bool operator!=(const Matrix& m) const
    {
        return !(*this == m);
    }

    // returns the number of rows
    int get_rows() const
    {
        return rows;
    }

    // returns the number of columns
    int get_cols() const
    {
        return cols;
    }

    // matrix(i, j) is the matrix's entry at the i'th row and j'th column, zero based (for get and set)
    T operator()(int i, int j) const
    {
        return major_row[i][j];
    }

    T& operator()(int i, int j)
    {
        return major_row[i][j];
    }

    // matrix at(i,j) returns the matrix's entry at the i'th row and j'th column, if there is such entry, zero based, or NULL if it does not exist.
    T at(int i, int j) const
    {
        if (i < rows && j < cols && i >= 0 && j >= 0) {
            return (*this)(i, j);
        }
        else {
            return NULL;
        }
    }

    // changes the matrix to its transposed form
    void transpose_in_place()
    {
        T temp;
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < i; j++) {
                temp = (*this)(i, j);
                (*this)(i, j) = (*this)(j, i);
                (*this)(j, i) = temp;
            }
        }
    }

    // returns the transposed matrix.
    Matrix transpose() const
    {
        Matrix transed(cols, rows);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                transed(j, i) = (*this)(i, j);
            }
        }
        return transed;
    }

    // in-place matrices addition
    void operator+=(const Matrix& m)
    {
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                (*this)(i, j) += m(i, j);
            }
        }
    }

    // matrices addition
    Matrix operator+(const Matrix& m) const
    {
        Matrix sum(rows, cols);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                sum(i, j) = (*this)(i, j) + m(i, j);
            }
        }
        return sum;
    }

    // matrices multiplication
    Matrix operator*(const Matrix& m) const
    {
        Matrix prod(rows, m.get_cols());
        T tmp_val;
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < m.get_cols(); ++j) {
                tmp_val = 0;
                for (int k = 0; k < cols; ++k) {
                    tmp_val += (*this)(i, k) * m(k, j);
                }
                prod(i, j) = tmp_val;
            }
        }
        return prod;
    }

    // two matrices with 2**power get two matrices A and B
    Matrix strassen(const Matrix<T>& B)
    {
        Matrix A = *this;
        if (A.get_cols() != B.get_rows()) {
            throw std::runtime_error("Very bad sizes!"); // todo
        }
        int power = round_up_power_2(std::max({A.get_rows(), A.get_cols(), B.get_cols()}));
        int paddedSize = 1 << power;

        // padding
        Matrix paddedA(paddedSize, paddedSize);
        Matrix paddedB(paddedSize, paddedSize);

        for (int i = 0; i < A.get_rows(); ++i) {
            for (int j = 0; j < A.get_cols(); ++j) {
                paddedA(i, j) = A(i, j);
            }
        }

        for (int i = 0; i < B.get_rows(); ++i) {
            for (int j = 0; j < B.get_cols(); ++j) {
                paddedB(i, j) = B(i, j);
            }
        }

        // calculate the product with mul()
        Matrix paddedProd = mul(paddedA, paddedB, power);

        // unpad the product
        Matrix prod(A.get_rows(), B.get_cols());
        for (int i = 0; i < prod.get_rows(); ++i) {
            for (int j = 0; j < prod.get_cols(); ++j) {
                prod(i, j) = paddedProd(i, j);
            }
        }

        return prod;
    }

    // matrix negation
    Matrix operator-() const
    {
        Matrix neg(rows, cols);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                neg(i, j) = -(*this)(i, j);
            }
        }
        return neg;
    }

    // in-place matrices substraction
    void operator-=(const Matrix& m)
    {
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; ++j) {
                (*this)(i, j) -= m(i, j);
            }
        }
    }

    // matrices substraction
    Matrix operator-(const Matrix& m) const
    {
        Matrix diff(rows, cols);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; ++j) {
                diff(i, j) = (*this)(i, j) - m(i, j);
            }
        }
        return diff;
    }

    // matrix by vector multiplication
    friend Vector<T> operator*(const Matrix& m, const Vector<T>& v)
    {
        Vector<T> res(m.get_rows());
        T tmp_val;
        for (int i = 0; i < m.get_rows(); ++i) {
            tmp_val = 0;
            for (int j = 0; j < m.get_cols(); ++j) {
                tmp_val += m(i, j) * v(j);
            }
            res(i) = tmp_val;
        }
        return res;
    }

    // in-place matrix by scalar multiplication
    void operator*=(const T& s)
    {
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                (*this)(i, j) *= s;
            }
        }
    }

    // matrix by scalar multiplication (both directions)
    friend Matrix operator*(const Matrix& m, const T& s)
    {
        Matrix prod(m.get_rows(), m.get_cols());
        for (int i = 0; i < m.get_rows(); ++i) {
            for (int j = 0; j < m.get_cols(); ++j) {
                prod(i, j) = m(i, j) * s;
            }
        }
        return prod;
    }

    friend Matrix operator*(const T& s, const Matrix& m)
    {
        return m * s;
    }
};

// sending to output stream using <<
template <class T>
std::ostream& operator<<(std::ostream& strm, const Matrix<T>& m)
{
    strm << "[";
    for (int i = 0; i < m.get_rows(); ++i) {
        strm << "(";
        for (int j = 0; j < m.get_cols() - 1; ++j) {
            strm << m(i, j) << ", ";
        }
        strm << m(i, m.get_cols() - 1) << ")\n";
    }
    strm << "]" << std::endl;
    return strm;
}

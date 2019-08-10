#pragma once

#include "Vector.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

template <class T>
class Matrix {
    using vector_t = std::vector<T>;

private:
    unsigned int rows; // num of rows
    unsigned int cols; // num of columns
    T determinant = 0;
    unsigned int rank = 0;
    bool det_updated = true;
    bool rank_updated = true;
    std::vector<vector_t> data; // the matrix's entries, stored as a vector of rows of the matrix, in row major fashion.

    // if matrix is square, the function computes the inverse if exists and updates the rank and det
    // if the matrix is not square, it only updates the rank
    Matrix gaussian_elimination()
    {
        rank = cols;
        determinant = 1;
        Matrix copy(*this);
        Matrix inverse(rows, rows);
        // set inverse to unit matrix
        for (unsigned int i = 0; i < rows; ++i) {
            inverse(i, i) = 1;
        }
        int first_non_zero = 0;
        for (unsigned int row = 0; row < rows && first_non_zero < cols;) {
            // if copy(i,i) is 0, swap it with a lower row.
            for (unsigned int lower_row = row + 1; lower_row < rows && copy(row, first_non_zero) == 0; ++lower_row) {
                if (copy(lower_row, first_non_zero) != 0) {
                    copy.swap_rows(row, lower_row);
                    inverse.swap_rows(row, lower_row);
                    determinant *= -1;
                }
            }

            if (copy(row, first_non_zero) != 0) {
                T scalar = copy(row, first_non_zero);
                copy.multiply_row_by_scalar(row, 1 / scalar);
                inverse.multiply_row_by_scalar(row, 1 / scalar);
                determinant *= scalar;
                // set all other rows to zero in the first_non_zero-th index
                for (unsigned int other_row = 0; other_row < rows; ++other_row) {
                    scalar = copy(other_row, first_non_zero);
                    if (other_row != row && scalar != 0) {
                        copy.add_multiplied_row(other_row, row, scalar * (-1));
                        inverse.add_multiplied_row(other_row, row, scalar * (-1));
                    }
                }
                ++first_non_zero;
                ++row;
            }
            else {
                determinant = 0;
                ++first_non_zero;
                --rank;
            }
        }

        rank_updated = true;
        det_updated = true;
        if (determinant != 0) {
            inverse.determinant = 1 / determinant;
            inverse.det_updated = true;
            inverse.rank = rows;
            inverse.rank_updated = true;
        }
        return std::move(inverse);
    }

    // max(s: 2^s<=x) x must be nonegative
    int round_up_power_2(int x)
    {
        int i = 0;
        while (x > (1 << i)) {
            ++i;
        }
        return i;
    }

    static void add_in_place(Matrix& A, const Matrix& B, int xA, int yA, int xB, int yB, int size)
    {
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                A(xA + i, yA + j) += B(xB + i, yB + j);
            }
        }
    }

    static void sub_in_place(Matrix& A, const Matrix& B, int xA, int yA, int xB, int yB, int size)
    {
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                A(xA + i, yA + j) -= B(xB + i, yB + j);
            }
        }
    }

    static void copy(Matrix& A, const Matrix& B, int xA, int yA, int xB, int yB, int size)
    {
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                A(xA + i, yA + j) = B(xB + i, yB + j);
            }
        }
    }

    static void base_case8(const Matrix& A, const Matrix& B, Matrix& C, int xA, int yA, int xB, int yB, int xC, int yC);

    static void base_case4(const Matrix& paddedA, const Matrix& paddedB, Matrix& paddedProd);

    static void mul(const Matrix& A, const Matrix& B, Matrix& C, int xA, int yA, int xB, int yB, int xC, int yC, int size)
    {
        if (size == 8) {
            Matrix::base_case8(A, B, C, xA, yA, xB, yB, xC, yC);
            return;
        }

        int half_size = size / 2;

        Matrix X(half_size, half_size);
        Matrix Y(half_size, half_size);
        Matrix M(half_size, half_size);

        // M1
        copy(X, A, 0, 0, xA, yA, half_size);
        add_in_place(X, A, 0, 0, xA + half_size, yA + half_size, half_size);
        copy(Y, B, 0, 0, xB, yB, half_size);
        add_in_place(Y, B, 0, 0, xB + half_size, yB + half_size, half_size);

        mul(X, Y, M, 0, 0, 0, 0, 0, 0, half_size);
        copy(C, M, xC + half_size, yC + half_size, 0, 0, half_size);
        copy(C, M, xC, yC, 0, 0, half_size);

        // M2
        copy(X, A, 0, 0, xA + half_size, yA, half_size);
        add_in_place(X, A, 0, 0, xA + half_size, yA + half_size, half_size);

        mul(X, B, M, 0, 0, xB, yB, 0, 0, half_size);
        copy(C, M, xC + half_size, yC, 0, 0, half_size);
        sub_in_place(C, M, xC + half_size, yC + half_size, 0, 0, half_size);

        // M3
        copy(Y, B, 0, 0, xB, yB + half_size, half_size);
        sub_in_place(Y, B, 0, 0, xB + half_size, yB + half_size, half_size);

        mul(A, Y, M, xA, yA, 0, 0, 0, 0, half_size);
        copy(C, M, xC, yC + half_size, 0, 0, half_size);
        add_in_place(C, M, xC + half_size, yC + half_size, 0, 0, half_size);

        // M4
        copy(Y, B, 0, 0, xB + half_size, yB, half_size);
        sub_in_place(Y, B, 0, 0, xB, yB, half_size);

        mul(A, Y, M, xA + half_size, yA + half_size, 0, 0, 0, 0, half_size);
        add_in_place(C, M, xC, yC, 0, 0, half_size);
        add_in_place(C, M, xC + half_size, yC, 0, 0, half_size);

        // M5
        copy(X, A, 0, 0, xA, yA, half_size);
        add_in_place(X, A, 0, 0, xA, yA + half_size, half_size);

        mul(X, B, M, 0, 0, xB + half_size, yB + half_size, 0, 0, half_size);
        sub_in_place(C, M, xC, yC, 0, 0, half_size);
        add_in_place(C, M, xC, yC + half_size, 0, 0, half_size);

        // M6
        copy(X, A, 0, 0, xA + half_size, yA, half_size);
        sub_in_place(X, A, 0, 0, xA, yA, half_size);
        copy(Y, B, 0, 0, xB, yB, half_size);
        add_in_place(Y, B, 0, 0, xB, yB + half_size, half_size);

        mul(X, Y, M, 0, 0, 0, 0, 0, 0, half_size);
        add_in_place(C, M, xC + half_size, yC + half_size, 0, 0, half_size);

        // M7
        copy(X, A, 0, 0, xA, yA + half_size, half_size);
        sub_in_place(X, A, 0, 0, xA + half_size, yA + half_size, half_size);
        copy(Y, B, 0, 0, xB + half_size, yB, half_size);
        add_in_place(Y, B, 0, 0, xB + half_size, yB + half_size, half_size);

        mul(X, Y, M, 0, 0, 0, 0, 0, 0, half_size);
        add_in_place(C, M, xC, yC, 0, 0, half_size);
    }

public:
    // default constructor. Constructs a 1X1 with 0
    explicit Matrix() : rows(1), cols(1), data(std::vector<vector_t>(1, vector_t(1, 0)))
    {
    }

    // constructs an rXc matrix, filled with 0's
    explicit Matrix(int r, int c) : rows(r), cols(c), data(std::vector<vector_t>(rows, vector_t(cols, 0)))
    {
    }

    // constructor by const reference to std::vector of vector_t (the entries)
    explicit Matrix(const std::vector<vector_t>& entries) : rows(entries.size()), cols(entries[0].size()), data(entries), det_updated(false), rank_updated(false)
    {
    }

    // constructor by rvalue reference to std::vector of vector_t (the entries)
    explicit Matrix(std::vector<vector_t>&& entries) : rows(entries.size()), cols(entries[0].size()), data(std::move(entries)), det_updated(false), rank_updated(false)
    {
    }

    // copy constructor. Constructs a copy of a given matrix
    Matrix(const Matrix& m) : rows(m.rows), cols(m.cols), data(m.data), det_updated(m.det_updated), rank_updated(m.rank_updated), determinant(m.determinant), rank(m.rank)
    {
    }

    // move constructor by rvalue reference of a different matrix
    Matrix(Matrix&& m) noexcept : rows(m.rows), cols(m.cols), data(std::move(m.data)), det_updated(m.det_updated), rank_updated(m.rank_updated), determinant(m.determinant), rank(m.rank)
    {
    }

    // copy assignment by const reference
    Matrix& operator=(const Matrix& m)
    {
        rows = m.rows;
        cols = m.cols;
        data = m.data;
        determinant = m.determinant;
        rank = m.rank;
        rank_updated = m.rank_updated;
        det_updated = m.det_updated;
    }

    // move assignment by rvalue reference // need to add col rows det and rank???
    Matrix& operator=(Matrix&& m) noexcept
    {
        data = std::move(m.data);
        rows = m.rows;
        cols = m.cols;
        determinant = m.determinant;
        rank = m.rank;
        rank_updated = m.rank_updated;
        det_updated = m.det_updated;
        return *this;
    }

    // return whether or not the given matrix has the same entries.
    bool operator==(const Matrix& m)
    {
        bool result = data == m.data;
        if (result) {
            if (det_updated && !m.det_updated) {
                m.determinant = determinant;
                m.det_updated = true;
            }
            else if (m.det_updated && !det_updated) {
                determinant = m.determinant;
                det_updated = true;
            }

            if (rank_updated && !m.rank_updated) {
                m.rank = rank;
                m.rank_updated = true;
            }
            else if (m.rank_updated && !rank_updated) {
                rank = m.rank;
                rank_updated = true;
            }
        }
        return result;
    }

    bool operator==(const Matrix& m) const
    {
        return data == m.data;
    }

    bool operator!=(const Matrix& m) const
    {
        return !(*this == m);
    }

    bool operator!=(const Matrix& m)
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
    const T& operator()(int i, int j) const
    {
        return data[i][j];
    }

    T& operator()(int i, int j)
    {
        rank_updated = false;
        det_updated = false;
        return data[i][j];
    }

    // matrix(i) is the matrix's i'th row, zero based (for get and set)
    const vector_t& operator()(unsigned int i) const
    {
        return data[i];
    }

    vector_t& operator()(unsigned int i)
    {
        det_updated = false;
        rank_updated = false;
        return data[i];
    }

    // changes the matrix to its transposed form
    // keeps the det and rank unchanged
    void transpose_in_place()
    {
        for (unsigned int i = 0; i < rows; ++i) {
            for (unsigned int j = 0; j < i; ++j) {
                std::swap((*this)(i, j), (*this)(j, i));
            }
        }
    }

    // returns the transposed matrix.
    // keeps the det and rank unchanged
    Matrix transpose() const
    {
        Matrix transed(cols, rows);
        for (unsigned int i = 0; i < rows; ++i) {
            for (unsigned int j = 0; j < cols; ++j) {
                transed(j, i) = (*this)(i, j);
            }
        }
        return std::move(transed);
    }

    // in-place matrices addition
    // det and rank are dirty
    Matrix& operator+=(const Matrix& m)
    {
        for (unsigned int i = 0; i < rows; ++i) {
            for (unsigned int j = 0; j < cols; ++j) {
                (*this)(i, j) += m(i, j);
            }
        }
        rank_updated = false;
        det_updated = false;
        return *this;
    }

    // matrices addition
    Matrix operator+(const Matrix& m) const
    {
        Matrix sum(rows, cols);
        for (unsigned int i = 0; i < rows; ++i) {
            for (unsigned int j = 0; j < cols; ++j) {
                sum(i, j) = (*this)(i, j) + m(i, j);
            }
        }
        return std::move(sum);
    }

    // matrices multiplication
    Matrix operator*(const Matrix& m) const
    {
        Matrix prod(rows, m.get_cols());
        for (unsigned int i = 0; i < rows; ++i) {
            for (unsigned int j = 0; j < m.cols; ++j) {
                prod(i, j) = 0;
                for (unsigned int k = 0; k < cols; ++k) {
                    prod(i, j) += (*this)(i, k) * m(k, j);
                }
            }
        }
        if (det_updated && m.det_updated) // can compute the new det
        {
            prod.determinant = m.determinant * determinant;
            prod.det_updated = true;
        }
        return std::move(prod);
    }

    // two matrices with 2**power get two matrices A and B
    Matrix strassen(const Matrix<T>& B)
    {
        Matrix A = *this;
        int power = round_up_power_2(std::max({A.get_rows(), A.get_cols(), B.get_cols()}));
        int size = 1 << power;

        // padding
        Matrix paddedA(size, size);
        Matrix paddedB(size, size);

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

        Matrix paddedProd = Matrix(size, size);
        if (size == 1) {
            paddedProd(0, 0) = A(0, 0) * B(0, 0);
        }

        else if (size == 2) {
            paddedProd(0, 0) = A(0, 0) * B(0, 0) + A(0, 1) * B(1, 0);
            paddedProd(0, 1) = A(0, 0) * B(0, 1) + A(0, 1) * B(1, 1);
            paddedProd(1, 0) = A(1, 0) * B(0, 0) + A(1, 1) * B(1, 0);
            paddedProd(1, 1) = A(1, 0) * B(0, 1) + A(1, 1) * B(1, 1);
        }

        else if (size == 4) {
            base_case4(paddedA, paddedB, paddedProd);
        }

        else {
            // calculate the product with mul()
            mul(paddedA, paddedB, paddedProd, 0, 0, 0, 0, 0, 0, size);
        }

        // unpad the product
        Matrix prod = Matrix(A.get_rows(), B.get_cols());
        for (int i = 0; i < prod.get_rows(); ++i) {
            for (int j = 0; j < prod.get_cols(); ++j) {
                prod(i, j) = paddedProd(i, j);
            }
        }

        if (det_updated && B.det_updated) // can compute the new det
        {
            prod.determinant = determinant * B.determinant;
            prod.det_updated = true;
        }

        return std::move(prod);
    }

    // matrix negation
    // rank unchanged, det can be updated accordingly
    Matrix operator-() const
    {
        return (*this) * (-1);
    }

    // in-place matrices substraction
    // rank and det unknown
    Matrix& operator-=(const Matrix& m)
    {
        for (unsigned int i = 0; i < rows; ++i) {
            for (unsigned int j = 0; j < cols; ++j) {
                (*this)(i, j) -= m(i, j);
            }
        }
        rank_updated = false;
        det_updated = false;
        return *this;
    }

    // matrices substraction
    Matrix operator-(const Matrix& m) const
    {
        Matrix diff(rows, cols);
        for (unsigned int i = 0; i < rows; ++i) {
            for (unsigned int j = 0; j < cols; ++j) {
                diff(i, j) = (*this)(i, j) - m(i, j);
            }
        }
        return std::move(diff);
    }

    // in-place matrix by scalar multiplication
    Matrix& operator*=(const T& s)
    {
        for (unsigned int i = 0; i < rows; ++i) {
            for (unsigned int j = 0; j < cols; ++j) {
                (*this)(i, j) *= s;
            }
        }
        if (det_updated) {
            determinant *= pow(s, rows);
        }
        if (s == 0) {
            rank = 0;
            rank_updated = true;
        }
        return *this;
    }

    // matrix by vector multiplication
    Vector<T> operator*(const Vector<T>& v) const
    {
        Vector<T> res(this->get_rows());
        T tmp_val;
        for (int i = 0; i < this->get_rows(); ++i) {
            tmp_val = 0;
            for (int j = 0; j < this->get_cols(); ++j) {
                tmp_val += (*this)(i, j) * v(j);
            }
            res(i) = tmp_val;
        }
        return std::move(res);
    }

    // matrix by scalar multiplication (both directions)
    friend Matrix operator*(const Matrix<T>& m, const T& s)
    {
        Matrix prod(m);
        prod *= s;
        return std::move(prod);
    }

    friend Matrix operator*(const T& s, const Matrix& m)
    {
        return std::move(m * s);
    }

    // gets two rows and swap them, the rows should have a valid range
    void swap_rows(int i, int j)
    {
        std::swap((*this)(i), (*this)(j));
        determinant *= -1;
    }

    void multiply_row_by_scalar(int row, T scalar)
    {
        for (unsigned int i = 0; i < cols; ++i) {
            (*this)(row, i) *= scalar;
        }
        determinant *= scalar;
        if (scalar == 0) {
            rank_updated = false;
        }
    }

    void add_multiplied_row(int row1, int row2, T scalar)
    {
        for (unsigned int i = 0; i < cols; ++i) {
            (*this)(row1, i) += (*this)(row2, i) * scalar;
        }
    }

    // gets one square matrix
    T get_det()
    {
        if (!det_updated) {
            gaussian_elimination();
        }
        return determinant;
    }
    // assumes the matrix has inverse
    Matrix find_inverse()
    {
        return gaussian_elimination();
    }

    unsigned int get_rank()
    {
        if (!rank_updated) {
            gaussian_elimination();
        }
        return rank;
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

template <class T>
inline void Matrix<T>::base_case8(const Matrix<T>& A, const Matrix<T>& B, Matrix<T>& C, int xA, int yA, int xB, int yB, int xC, int yC)
{
    C(xC + 0, yC + 0) = A(xA + 0, yA + 0) * B(xB + 0, yB + 0) + A(xA + 0, yA + 1) * B(xB + 1, yB + 0) + A(xA + 0, yA + 2) * B(xB + 2, yB + 0) + A(xA + 0, yA + 3) * B(xB + 3, yB + 0) + A(xA + 0, yA + 4) * B(xB + 4, yB + 0) + A(xA + 0, yA + 5) * B(xB + 5, yB + 0) + A(xA + 0, yA + 6) * B(xB + 6, yB + 0) + A(xA + 0, yA + 7) * B(xB + 7, yB + 0);
    C(xC + 0, yC + 1) = A(xA + 0, yA + 0) * B(xB + 0, yB + 1) + A(xA + 0, yA + 1) * B(xB + 1, yB + 1) + A(xA + 0, yA + 2) * B(xB + 2, yB + 1) + A(xA + 0, yA + 3) * B(xB + 3, yB + 1) + A(xA + 0, yA + 4) * B(xB + 4, yB + 1) + A(xA + 0, yA + 5) * B(xB + 5, yB + 1) + A(xA + 0, yA + 6) * B(xB + 6, yB + 1) + A(xA + 0, yA + 7) * B(xB + 7, yB + 1);
    C(xC + 0, yC + 2) = A(xA + 0, yA + 0) * B(xB + 0, yB + 2) + A(xA + 0, yA + 1) * B(xB + 1, yB + 2) + A(xA + 0, yA + 2) * B(xB + 2, yB + 2) + A(xA + 0, yA + 3) * B(xB + 3, yB + 2) + A(xA + 0, yA + 4) * B(xB + 4, yB + 2) + A(xA + 0, yA + 5) * B(xB + 5, yB + 2) + A(xA + 0, yA + 6) * B(xB + 6, yB + 2) + A(xA + 0, yA + 7) * B(xB + 7, yB + 2);
    C(xC + 0, yC + 3) = A(xA + 0, yA + 0) * B(xB + 0, yB + 3) + A(xA + 0, yA + 1) * B(xB + 1, yB + 3) + A(xA + 0, yA + 2) * B(xB + 2, yB + 3) + A(xA + 0, yA + 3) * B(xB + 3, yB + 3) + A(xA + 0, yA + 4) * B(xB + 4, yB + 3) + A(xA + 0, yA + 5) * B(xB + 5, yB + 3) + A(xA + 0, yA + 6) * B(xB + 6, yB + 3) + A(xA + 0, yA + 7) * B(xB + 7, yB + 3);
    C(xC + 0, yC + 4) = A(xA + 0, yA + 0) * B(xB + 0, yB + 4) + A(xA + 0, yA + 1) * B(xB + 1, yB + 4) + A(xA + 0, yA + 2) * B(xB + 2, yB + 4) + A(xA + 0, yA + 3) * B(xB + 3, yB + 4) + A(xA + 0, yA + 4) * B(xB + 4, yB + 4) + A(xA + 0, yA + 5) * B(xB + 5, yB + 4) + A(xA + 0, yA + 6) * B(xB + 6, yB + 4) + A(xA + 0, yA + 7) * B(xB + 7, yB + 4);
    C(xC + 0, yC + 5) = A(xA + 0, yA + 0) * B(xB + 0, yB + 5) + A(xA + 0, yA + 1) * B(xB + 1, yB + 5) + A(xA + 0, yA + 2) * B(xB + 2, yB + 5) + A(xA + 0, yA + 3) * B(xB + 3, yB + 5) + A(xA + 0, yA + 4) * B(xB + 4, yB + 5) + A(xA + 0, yA + 5) * B(xB + 5, yB + 5) + A(xA + 0, yA + 6) * B(xB + 6, yB + 5) + A(xA + 0, yA + 7) * B(xB + 7, yB + 5);
    C(xC + 0, yC + 6) = A(xA + 0, yA + 0) * B(xB + 0, yB + 6) + A(xA + 0, yA + 1) * B(xB + 1, yB + 6) + A(xA + 0, yA + 2) * B(xB + 2, yB + 6) + A(xA + 0, yA + 3) * B(xB + 3, yB + 6) + A(xA + 0, yA + 4) * B(xB + 4, yB + 6) + A(xA + 0, yA + 5) * B(xB + 5, yB + 6) + A(xA + 0, yA + 6) * B(xB + 6, yB + 6) + A(xA + 0, yA + 7) * B(xB + 7, yB + 6);
    C(xC + 0, yC + 7) = A(xA + 0, yA + 0) * B(xB + 0, yB + 7) + A(xA + 0, yA + 1) * B(xB + 1, yB + 7) + A(xA + 0, yA + 2) * B(xB + 2, yB + 7) + A(xA + 0, yA + 3) * B(xB + 3, yB + 7) + A(xA + 0, yA + 4) * B(xB + 4, yB + 7) + A(xA + 0, yA + 5) * B(xB + 5, yB + 7) + A(xA + 0, yA + 6) * B(xB + 6, yB + 7) + A(xA + 0, yA + 7) * B(xB + 7, yB + 7);
    C(xC + 1, yC + 0) = A(xA + 1, yA + 0) * B(xB + 0, yB + 0) + A(xA + 1, yA + 1) * B(xB + 1, yB + 0) + A(xA + 1, yA + 2) * B(xB + 2, yB + 0) + A(xA + 1, yA + 3) * B(xB + 3, yB + 0) + A(xA + 1, yA + 4) * B(xB + 4, yB + 0) + A(xA + 1, yA + 5) * B(xB + 5, yB + 0) + A(xA + 1, yA + 6) * B(xB + 6, yB + 0) + A(xA + 1, yA + 7) * B(xB + 7, yB + 0);
    C(xC + 1, yC + 1) = A(xA + 1, yA + 0) * B(xB + 0, yB + 1) + A(xA + 1, yA + 1) * B(xB + 1, yB + 1) + A(xA + 1, yA + 2) * B(xB + 2, yB + 1) + A(xA + 1, yA + 3) * B(xB + 3, yB + 1) + A(xA + 1, yA + 4) * B(xB + 4, yB + 1) + A(xA + 1, yA + 5) * B(xB + 5, yB + 1) + A(xA + 1, yA + 6) * B(xB + 6, yB + 1) + A(xA + 1, yA + 7) * B(xB + 7, yB + 1);
    C(xC + 1, yC + 2) = A(xA + 1, yA + 0) * B(xB + 0, yB + 2) + A(xA + 1, yA + 1) * B(xB + 1, yB + 2) + A(xA + 1, yA + 2) * B(xB + 2, yB + 2) + A(xA + 1, yA + 3) * B(xB + 3, yB + 2) + A(xA + 1, yA + 4) * B(xB + 4, yB + 2) + A(xA + 1, yA + 5) * B(xB + 5, yB + 2) + A(xA + 1, yA + 6) * B(xB + 6, yB + 2) + A(xA + 1, yA + 7) * B(xB + 7, yB + 2);
    C(xC + 1, yC + 3) = A(xA + 1, yA + 0) * B(xB + 0, yB + 3) + A(xA + 1, yA + 1) * B(xB + 1, yB + 3) + A(xA + 1, yA + 2) * B(xB + 2, yB + 3) + A(xA + 1, yA + 3) * B(xB + 3, yB + 3) + A(xA + 1, yA + 4) * B(xB + 4, yB + 3) + A(xA + 1, yA + 5) * B(xB + 5, yB + 3) + A(xA + 1, yA + 6) * B(xB + 6, yB + 3) + A(xA + 1, yA + 7) * B(xB + 7, yB + 3);
    C(xC + 1, yC + 4) = A(xA + 1, yA + 0) * B(xB + 0, yB + 4) + A(xA + 1, yA + 1) * B(xB + 1, yB + 4) + A(xA + 1, yA + 2) * B(xB + 2, yB + 4) + A(xA + 1, yA + 3) * B(xB + 3, yB + 4) + A(xA + 1, yA + 4) * B(xB + 4, yB + 4) + A(xA + 1, yA + 5) * B(xB + 5, yB + 4) + A(xA + 1, yA + 6) * B(xB + 6, yB + 4) + A(xA + 1, yA + 7) * B(xB + 7, yB + 4);
    C(xC + 1, yC + 5) = A(xA + 1, yA + 0) * B(xB + 0, yB + 5) + A(xA + 1, yA + 1) * B(xB + 1, yB + 5) + A(xA + 1, yA + 2) * B(xB + 2, yB + 5) + A(xA + 1, yA + 3) * B(xB + 3, yB + 5) + A(xA + 1, yA + 4) * B(xB + 4, yB + 5) + A(xA + 1, yA + 5) * B(xB + 5, yB + 5) + A(xA + 1, yA + 6) * B(xB + 6, yB + 5) + A(xA + 1, yA + 7) * B(xB + 7, yB + 5);
    C(xC + 1, yC + 6) = A(xA + 1, yA + 0) * B(xB + 0, yB + 6) + A(xA + 1, yA + 1) * B(xB + 1, yB + 6) + A(xA + 1, yA + 2) * B(xB + 2, yB + 6) + A(xA + 1, yA + 3) * B(xB + 3, yB + 6) + A(xA + 1, yA + 4) * B(xB + 4, yB + 6) + A(xA + 1, yA + 5) * B(xB + 5, yB + 6) + A(xA + 1, yA + 6) * B(xB + 6, yB + 6) + A(xA + 1, yA + 7) * B(xB + 7, yB + 6);
    C(xC + 1, yC + 7) = A(xA + 1, yA + 0) * B(xB + 0, yB + 7) + A(xA + 1, yA + 1) * B(xB + 1, yB + 7) + A(xA + 1, yA + 2) * B(xB + 2, yB + 7) + A(xA + 1, yA + 3) * B(xB + 3, yB + 7) + A(xA + 1, yA + 4) * B(xB + 4, yB + 7) + A(xA + 1, yA + 5) * B(xB + 5, yB + 7) + A(xA + 1, yA + 6) * B(xB + 6, yB + 7) + A(xA + 1, yA + 7) * B(xB + 7, yB + 7);
    C(xC + 2, yC + 0) = A(xA + 2, yA + 0) * B(xB + 0, yB + 0) + A(xA + 2, yA + 1) * B(xB + 1, yB + 0) + A(xA + 2, yA + 2) * B(xB + 2, yB + 0) + A(xA + 2, yA + 3) * B(xB + 3, yB + 0) + A(xA + 2, yA + 4) * B(xB + 4, yB + 0) + A(xA + 2, yA + 5) * B(xB + 5, yB + 0) + A(xA + 2, yA + 6) * B(xB + 6, yB + 0) + A(xA + 2, yA + 7) * B(xB + 7, yB + 0);
    C(xC + 2, yC + 1) = A(xA + 2, yA + 0) * B(xB + 0, yB + 1) + A(xA + 2, yA + 1) * B(xB + 1, yB + 1) + A(xA + 2, yA + 2) * B(xB + 2, yB + 1) + A(xA + 2, yA + 3) * B(xB + 3, yB + 1) + A(xA + 2, yA + 4) * B(xB + 4, yB + 1) + A(xA + 2, yA + 5) * B(xB + 5, yB + 1) + A(xA + 2, yA + 6) * B(xB + 6, yB + 1) + A(xA + 2, yA + 7) * B(xB + 7, yB + 1);
    C(xC + 2, yC + 2) = A(xA + 2, yA + 0) * B(xB + 0, yB + 2) + A(xA + 2, yA + 1) * B(xB + 1, yB + 2) + A(xA + 2, yA + 2) * B(xB + 2, yB + 2) + A(xA + 2, yA + 3) * B(xB + 3, yB + 2) + A(xA + 2, yA + 4) * B(xB + 4, yB + 2) + A(xA + 2, yA + 5) * B(xB + 5, yB + 2) + A(xA + 2, yA + 6) * B(xB + 6, yB + 2) + A(xA + 2, yA + 7) * B(xB + 7, yB + 2);
    C(xC + 2, yC + 3) = A(xA + 2, yA + 0) * B(xB + 0, yB + 3) + A(xA + 2, yA + 1) * B(xB + 1, yB + 3) + A(xA + 2, yA + 2) * B(xB + 2, yB + 3) + A(xA + 2, yA + 3) * B(xB + 3, yB + 3) + A(xA + 2, yA + 4) * B(xB + 4, yB + 3) + A(xA + 2, yA + 5) * B(xB + 5, yB + 3) + A(xA + 2, yA + 6) * B(xB + 6, yB + 3) + A(xA + 2, yA + 7) * B(xB + 7, yB + 3);
    C(xC + 2, yC + 4) = A(xA + 2, yA + 0) * B(xB + 0, yB + 4) + A(xA + 2, yA + 1) * B(xB + 1, yB + 4) + A(xA + 2, yA + 2) * B(xB + 2, yB + 4) + A(xA + 2, yA + 3) * B(xB + 3, yB + 4) + A(xA + 2, yA + 4) * B(xB + 4, yB + 4) + A(xA + 2, yA + 5) * B(xB + 5, yB + 4) + A(xA + 2, yA + 6) * B(xB + 6, yB + 4) + A(xA + 2, yA + 7) * B(xB + 7, yB + 4);
    C(xC + 2, yC + 5) = A(xA + 2, yA + 0) * B(xB + 0, yB + 5) + A(xA + 2, yA + 1) * B(xB + 1, yB + 5) + A(xA + 2, yA + 2) * B(xB + 2, yB + 5) + A(xA + 2, yA + 3) * B(xB + 3, yB + 5) + A(xA + 2, yA + 4) * B(xB + 4, yB + 5) + A(xA + 2, yA + 5) * B(xB + 5, yB + 5) + A(xA + 2, yA + 6) * B(xB + 6, yB + 5) + A(xA + 2, yA + 7) * B(xB + 7, yB + 5);
    C(xC + 2, yC + 6) = A(xA + 2, yA + 0) * B(xB + 0, yB + 6) + A(xA + 2, yA + 1) * B(xB + 1, yB + 6) + A(xA + 2, yA + 2) * B(xB + 2, yB + 6) + A(xA + 2, yA + 3) * B(xB + 3, yB + 6) + A(xA + 2, yA + 4) * B(xB + 4, yB + 6) + A(xA + 2, yA + 5) * B(xB + 5, yB + 6) + A(xA + 2, yA + 6) * B(xB + 6, yB + 6) + A(xA + 2, yA + 7) * B(xB + 7, yB + 6);
    C(xC + 2, yC + 7) = A(xA + 2, yA + 0) * B(xB + 0, yB + 7) + A(xA + 2, yA + 1) * B(xB + 1, yB + 7) + A(xA + 2, yA + 2) * B(xB + 2, yB + 7) + A(xA + 2, yA + 3) * B(xB + 3, yB + 7) + A(xA + 2, yA + 4) * B(xB + 4, yB + 7) + A(xA + 2, yA + 5) * B(xB + 5, yB + 7) + A(xA + 2, yA + 6) * B(xB + 6, yB + 7) + A(xA + 2, yA + 7) * B(xB + 7, yB + 7);
    C(xC + 3, yC + 0) = A(xA + 3, yA + 0) * B(xB + 0, yB + 0) + A(xA + 3, yA + 1) * B(xB + 1, yB + 0) + A(xA + 3, yA + 2) * B(xB + 2, yB + 0) + A(xA + 3, yA + 3) * B(xB + 3, yB + 0) + A(xA + 3, yA + 4) * B(xB + 4, yB + 0) + A(xA + 3, yA + 5) * B(xB + 5, yB + 0) + A(xA + 3, yA + 6) * B(xB + 6, yB + 0) + A(xA + 3, yA + 7) * B(xB + 7, yB + 0);
    C(xC + 3, yC + 1) = A(xA + 3, yA + 0) * B(xB + 0, yB + 1) + A(xA + 3, yA + 1) * B(xB + 1, yB + 1) + A(xA + 3, yA + 2) * B(xB + 2, yB + 1) + A(xA + 3, yA + 3) * B(xB + 3, yB + 1) + A(xA + 3, yA + 4) * B(xB + 4, yB + 1) + A(xA + 3, yA + 5) * B(xB + 5, yB + 1) + A(xA + 3, yA + 6) * B(xB + 6, yB + 1) + A(xA + 3, yA + 7) * B(xB + 7, yB + 1);
    C(xC + 3, yC + 2) = A(xA + 3, yA + 0) * B(xB + 0, yB + 2) + A(xA + 3, yA + 1) * B(xB + 1, yB + 2) + A(xA + 3, yA + 2) * B(xB + 2, yB + 2) + A(xA + 3, yA + 3) * B(xB + 3, yB + 2) + A(xA + 3, yA + 4) * B(xB + 4, yB + 2) + A(xA + 3, yA + 5) * B(xB + 5, yB + 2) + A(xA + 3, yA + 6) * B(xB + 6, yB + 2) + A(xA + 3, yA + 7) * B(xB + 7, yB + 2);
    C(xC + 3, yC + 3) = A(xA + 3, yA + 0) * B(xB + 0, yB + 3) + A(xA + 3, yA + 1) * B(xB + 1, yB + 3) + A(xA + 3, yA + 2) * B(xB + 2, yB + 3) + A(xA + 3, yA + 3) * B(xB + 3, yB + 3) + A(xA + 3, yA + 4) * B(xB + 4, yB + 3) + A(xA + 3, yA + 5) * B(xB + 5, yB + 3) + A(xA + 3, yA + 6) * B(xB + 6, yB + 3) + A(xA + 3, yA + 7) * B(xB + 7, yB + 3);
    C(xC + 3, yC + 4) = A(xA + 3, yA + 0) * B(xB + 0, yB + 4) + A(xA + 3, yA + 1) * B(xB + 1, yB + 4) + A(xA + 3, yA + 2) * B(xB + 2, yB + 4) + A(xA + 3, yA + 3) * B(xB + 3, yB + 4) + A(xA + 3, yA + 4) * B(xB + 4, yB + 4) + A(xA + 3, yA + 5) * B(xB + 5, yB + 4) + A(xA + 3, yA + 6) * B(xB + 6, yB + 4) + A(xA + 3, yA + 7) * B(xB + 7, yB + 4);
    C(xC + 3, yC + 5) = A(xA + 3, yA + 0) * B(xB + 0, yB + 5) + A(xA + 3, yA + 1) * B(xB + 1, yB + 5) + A(xA + 3, yA + 2) * B(xB + 2, yB + 5) + A(xA + 3, yA + 3) * B(xB + 3, yB + 5) + A(xA + 3, yA + 4) * B(xB + 4, yB + 5) + A(xA + 3, yA + 5) * B(xB + 5, yB + 5) + A(xA + 3, yA + 6) * B(xB + 6, yB + 5) + A(xA + 3, yA + 7) * B(xB + 7, yB + 5);
    C(xC + 3, yC + 6) = A(xA + 3, yA + 0) * B(xB + 0, yB + 6) + A(xA + 3, yA + 1) * B(xB + 1, yB + 6) + A(xA + 3, yA + 2) * B(xB + 2, yB + 6) + A(xA + 3, yA + 3) * B(xB + 3, yB + 6) + A(xA + 3, yA + 4) * B(xB + 4, yB + 6) + A(xA + 3, yA + 5) * B(xB + 5, yB + 6) + A(xA + 3, yA + 6) * B(xB + 6, yB + 6) + A(xA + 3, yA + 7) * B(xB + 7, yB + 6);
    C(xC + 3, yC + 7) = A(xA + 3, yA + 0) * B(xB + 0, yB + 7) + A(xA + 3, yA + 1) * B(xB + 1, yB + 7) + A(xA + 3, yA + 2) * B(xB + 2, yB + 7) + A(xA + 3, yA + 3) * B(xB + 3, yB + 7) + A(xA + 3, yA + 4) * B(xB + 4, yB + 7) + A(xA + 3, yA + 5) * B(xB + 5, yB + 7) + A(xA + 3, yA + 6) * B(xB + 6, yB + 7) + A(xA + 3, yA + 7) * B(xB + 7, yB + 7);
    C(xC + 4, yC + 0) = A(xA + 4, yA + 0) * B(xB + 0, yB + 0) + A(xA + 4, yA + 1) * B(xB + 1, yB + 0) + A(xA + 4, yA + 2) * B(xB + 2, yB + 0) + A(xA + 4, yA + 3) * B(xB + 3, yB + 0) + A(xA + 4, yA + 4) * B(xB + 4, yB + 0) + A(xA + 4, yA + 5) * B(xB + 5, yB + 0) + A(xA + 4, yA + 6) * B(xB + 6, yB + 0) + A(xA + 4, yA + 7) * B(xB + 7, yB + 0);
    C(xC + 4, yC + 1) = A(xA + 4, yA + 0) * B(xB + 0, yB + 1) + A(xA + 4, yA + 1) * B(xB + 1, yB + 1) + A(xA + 4, yA + 2) * B(xB + 2, yB + 1) + A(xA + 4, yA + 3) * B(xB + 3, yB + 1) + A(xA + 4, yA + 4) * B(xB + 4, yB + 1) + A(xA + 4, yA + 5) * B(xB + 5, yB + 1) + A(xA + 4, yA + 6) * B(xB + 6, yB + 1) + A(xA + 4, yA + 7) * B(xB + 7, yB + 1);
    C(xC + 4, yC + 2) = A(xA + 4, yA + 0) * B(xB + 0, yB + 2) + A(xA + 4, yA + 1) * B(xB + 1, yB + 2) + A(xA + 4, yA + 2) * B(xB + 2, yB + 2) + A(xA + 4, yA + 3) * B(xB + 3, yB + 2) + A(xA + 4, yA + 4) * B(xB + 4, yB + 2) + A(xA + 4, yA + 5) * B(xB + 5, yB + 2) + A(xA + 4, yA + 6) * B(xB + 6, yB + 2) + A(xA + 4, yA + 7) * B(xB + 7, yB + 2);
    C(xC + 4, yC + 3) = A(xA + 4, yA + 0) * B(xB + 0, yB + 3) + A(xA + 4, yA + 1) * B(xB + 1, yB + 3) + A(xA + 4, yA + 2) * B(xB + 2, yB + 3) + A(xA + 4, yA + 3) * B(xB + 3, yB + 3) + A(xA + 4, yA + 4) * B(xB + 4, yB + 3) + A(xA + 4, yA + 5) * B(xB + 5, yB + 3) + A(xA + 4, yA + 6) * B(xB + 6, yB + 3) + A(xA + 4, yA + 7) * B(xB + 7, yB + 3);
    C(xC + 4, yC + 4) = A(xA + 4, yA + 0) * B(xB + 0, yB + 4) + A(xA + 4, yA + 1) * B(xB + 1, yB + 4) + A(xA + 4, yA + 2) * B(xB + 2, yB + 4) + A(xA + 4, yA + 3) * B(xB + 3, yB + 4) + A(xA + 4, yA + 4) * B(xB + 4, yB + 4) + A(xA + 4, yA + 5) * B(xB + 5, yB + 4) + A(xA + 4, yA + 6) * B(xB + 6, yB + 4) + A(xA + 4, yA + 7) * B(xB + 7, yB + 4);
    C(xC + 4, yC + 5) = A(xA + 4, yA + 0) * B(xB + 0, yB + 5) + A(xA + 4, yA + 1) * B(xB + 1, yB + 5) + A(xA + 4, yA + 2) * B(xB + 2, yB + 5) + A(xA + 4, yA + 3) * B(xB + 3, yB + 5) + A(xA + 4, yA + 4) * B(xB + 4, yB + 5) + A(xA + 4, yA + 5) * B(xB + 5, yB + 5) + A(xA + 4, yA + 6) * B(xB + 6, yB + 5) + A(xA + 4, yA + 7) * B(xB + 7, yB + 5);
    C(xC + 4, yC + 6) = A(xA + 4, yA + 0) * B(xB + 0, yB + 6) + A(xA + 4, yA + 1) * B(xB + 1, yB + 6) + A(xA + 4, yA + 2) * B(xB + 2, yB + 6) + A(xA + 4, yA + 3) * B(xB + 3, yB + 6) + A(xA + 4, yA + 4) * B(xB + 4, yB + 6) + A(xA + 4, yA + 5) * B(xB + 5, yB + 6) + A(xA + 4, yA + 6) * B(xB + 6, yB + 6) + A(xA + 4, yA + 7) * B(xB + 7, yB + 6);
    C(xC + 4, yC + 7) = A(xA + 4, yA + 0) * B(xB + 0, yB + 7) + A(xA + 4, yA + 1) * B(xB + 1, yB + 7) + A(xA + 4, yA + 2) * B(xB + 2, yB + 7) + A(xA + 4, yA + 3) * B(xB + 3, yB + 7) + A(xA + 4, yA + 4) * B(xB + 4, yB + 7) + A(xA + 4, yA + 5) * B(xB + 5, yB + 7) + A(xA + 4, yA + 6) * B(xB + 6, yB + 7) + A(xA + 4, yA + 7) * B(xB + 7, yB + 7);
    C(xC + 5, yC + 0) = A(xA + 5, yA + 0) * B(xB + 0, yB + 0) + A(xA + 5, yA + 1) * B(xB + 1, yB + 0) + A(xA + 5, yA + 2) * B(xB + 2, yB + 0) + A(xA + 5, yA + 3) * B(xB + 3, yB + 0) + A(xA + 5, yA + 4) * B(xB + 4, yB + 0) + A(xA + 5, yA + 5) * B(xB + 5, yB + 0) + A(xA + 5, yA + 6) * B(xB + 6, yB + 0) + A(xA + 5, yA + 7) * B(xB + 7, yB + 0);
    C(xC + 5, yC + 1) = A(xA + 5, yA + 0) * B(xB + 0, yB + 1) + A(xA + 5, yA + 1) * B(xB + 1, yB + 1) + A(xA + 5, yA + 2) * B(xB + 2, yB + 1) + A(xA + 5, yA + 3) * B(xB + 3, yB + 1) + A(xA + 5, yA + 4) * B(xB + 4, yB + 1) + A(xA + 5, yA + 5) * B(xB + 5, yB + 1) + A(xA + 5, yA + 6) * B(xB + 6, yB + 1) + A(xA + 5, yA + 7) * B(xB + 7, yB + 1);
    C(xC + 5, yC + 2) = A(xA + 5, yA + 0) * B(xB + 0, yB + 2) + A(xA + 5, yA + 1) * B(xB + 1, yB + 2) + A(xA + 5, yA + 2) * B(xB + 2, yB + 2) + A(xA + 5, yA + 3) * B(xB + 3, yB + 2) + A(xA + 5, yA + 4) * B(xB + 4, yB + 2) + A(xA + 5, yA + 5) * B(xB + 5, yB + 2) + A(xA + 5, yA + 6) * B(xB + 6, yB + 2) + A(xA + 5, yA + 7) * B(xB + 7, yB + 2);
    C(xC + 5, yC + 3) = A(xA + 5, yA + 0) * B(xB + 0, yB + 3) + A(xA + 5, yA + 1) * B(xB + 1, yB + 3) + A(xA + 5, yA + 2) * B(xB + 2, yB + 3) + A(xA + 5, yA + 3) * B(xB + 3, yB + 3) + A(xA + 5, yA + 4) * B(xB + 4, yB + 3) + A(xA + 5, yA + 5) * B(xB + 5, yB + 3) + A(xA + 5, yA + 6) * B(xB + 6, yB + 3) + A(xA + 5, yA + 7) * B(xB + 7, yB + 3);
    C(xC + 5, yC + 4) = A(xA + 5, yA + 0) * B(xB + 0, yB + 4) + A(xA + 5, yA + 1) * B(xB + 1, yB + 4) + A(xA + 5, yA + 2) * B(xB + 2, yB + 4) + A(xA + 5, yA + 3) * B(xB + 3, yB + 4) + A(xA + 5, yA + 4) * B(xB + 4, yB + 4) + A(xA + 5, yA + 5) * B(xB + 5, yB + 4) + A(xA + 5, yA + 6) * B(xB + 6, yB + 4) + A(xA + 5, yA + 7) * B(xB + 7, yB + 4);
    C(xC + 5, yC + 5) = A(xA + 5, yA + 0) * B(xB + 0, yB + 5) + A(xA + 5, yA + 1) * B(xB + 1, yB + 5) + A(xA + 5, yA + 2) * B(xB + 2, yB + 5) + A(xA + 5, yA + 3) * B(xB + 3, yB + 5) + A(xA + 5, yA + 4) * B(xB + 4, yB + 5) + A(xA + 5, yA + 5) * B(xB + 5, yB + 5) + A(xA + 5, yA + 6) * B(xB + 6, yB + 5) + A(xA + 5, yA + 7) * B(xB + 7, yB + 5);
    C(xC + 5, yC + 6) = A(xA + 5, yA + 0) * B(xB + 0, yB + 6) + A(xA + 5, yA + 1) * B(xB + 1, yB + 6) + A(xA + 5, yA + 2) * B(xB + 2, yB + 6) + A(xA + 5, yA + 3) * B(xB + 3, yB + 6) + A(xA + 5, yA + 4) * B(xB + 4, yB + 6) + A(xA + 5, yA + 5) * B(xB + 5, yB + 6) + A(xA + 5, yA + 6) * B(xB + 6, yB + 6) + A(xA + 5, yA + 7) * B(xB + 7, yB + 6);
    C(xC + 5, yC + 7) = A(xA + 5, yA + 0) * B(xB + 0, yB + 7) + A(xA + 5, yA + 1) * B(xB + 1, yB + 7) + A(xA + 5, yA + 2) * B(xB + 2, yB + 7) + A(xA + 5, yA + 3) * B(xB + 3, yB + 7) + A(xA + 5, yA + 4) * B(xB + 4, yB + 7) + A(xA + 5, yA + 5) * B(xB + 5, yB + 7) + A(xA + 5, yA + 6) * B(xB + 6, yB + 7) + A(xA + 5, yA + 7) * B(xB + 7, yB + 7);
    C(xC + 6, yC + 0) = A(xA + 6, yA + 0) * B(xB + 0, yB + 0) + A(xA + 6, yA + 1) * B(xB + 1, yB + 0) + A(xA + 6, yA + 2) * B(xB + 2, yB + 0) + A(xA + 6, yA + 3) * B(xB + 3, yB + 0) + A(xA + 6, yA + 4) * B(xB + 4, yB + 0) + A(xA + 6, yA + 5) * B(xB + 5, yB + 0) + A(xA + 6, yA + 6) * B(xB + 6, yB + 0) + A(xA + 6, yA + 7) * B(xB + 7, yB + 0);
    C(xC + 6, yC + 1) = A(xA + 6, yA + 0) * B(xB + 0, yB + 1) + A(xA + 6, yA + 1) * B(xB + 1, yB + 1) + A(xA + 6, yA + 2) * B(xB + 2, yB + 1) + A(xA + 6, yA + 3) * B(xB + 3, yB + 1) + A(xA + 6, yA + 4) * B(xB + 4, yB + 1) + A(xA + 6, yA + 5) * B(xB + 5, yB + 1) + A(xA + 6, yA + 6) * B(xB + 6, yB + 1) + A(xA + 6, yA + 7) * B(xB + 7, yB + 1);
    C(xC + 6, yC + 2) = A(xA + 6, yA + 0) * B(xB + 0, yB + 2) + A(xA + 6, yA + 1) * B(xB + 1, yB + 2) + A(xA + 6, yA + 2) * B(xB + 2, yB + 2) + A(xA + 6, yA + 3) * B(xB + 3, yB + 2) + A(xA + 6, yA + 4) * B(xB + 4, yB + 2) + A(xA + 6, yA + 5) * B(xB + 5, yB + 2) + A(xA + 6, yA + 6) * B(xB + 6, yB + 2) + A(xA + 6, yA + 7) * B(xB + 7, yB + 2);
    C(xC + 6, yC + 3) = A(xA + 6, yA + 0) * B(xB + 0, yB + 3) + A(xA + 6, yA + 1) * B(xB + 1, yB + 3) + A(xA + 6, yA + 2) * B(xB + 2, yB + 3) + A(xA + 6, yA + 3) * B(xB + 3, yB + 3) + A(xA + 6, yA + 4) * B(xB + 4, yB + 3) + A(xA + 6, yA + 5) * B(xB + 5, yB + 3) + A(xA + 6, yA + 6) * B(xB + 6, yB + 3) + A(xA + 6, yA + 7) * B(xB + 7, yB + 3);
    C(xC + 6, yC + 4) = A(xA + 6, yA + 0) * B(xB + 0, yB + 4) + A(xA + 6, yA + 1) * B(xB + 1, yB + 4) + A(xA + 6, yA + 2) * B(xB + 2, yB + 4) + A(xA + 6, yA + 3) * B(xB + 3, yB + 4) + A(xA + 6, yA + 4) * B(xB + 4, yB + 4) + A(xA + 6, yA + 5) * B(xB + 5, yB + 4) + A(xA + 6, yA + 6) * B(xB + 6, yB + 4) + A(xA + 6, yA + 7) * B(xB + 7, yB + 4);
    C(xC + 6, yC + 5) = A(xA + 6, yA + 0) * B(xB + 0, yB + 5) + A(xA + 6, yA + 1) * B(xB + 1, yB + 5) + A(xA + 6, yA + 2) * B(xB + 2, yB + 5) + A(xA + 6, yA + 3) * B(xB + 3, yB + 5) + A(xA + 6, yA + 4) * B(xB + 4, yB + 5) + A(xA + 6, yA + 5) * B(xB + 5, yB + 5) + A(xA + 6, yA + 6) * B(xB + 6, yB + 5) + A(xA + 6, yA + 7) * B(xB + 7, yB + 5);
    C(xC + 6, yC + 6) = A(xA + 6, yA + 0) * B(xB + 0, yB + 6) + A(xA + 6, yA + 1) * B(xB + 1, yB + 6) + A(xA + 6, yA + 2) * B(xB + 2, yB + 6) + A(xA + 6, yA + 3) * B(xB + 3, yB + 6) + A(xA + 6, yA + 4) * B(xB + 4, yB + 6) + A(xA + 6, yA + 5) * B(xB + 5, yB + 6) + A(xA + 6, yA + 6) * B(xB + 6, yB + 6) + A(xA + 6, yA + 7) * B(xB + 7, yB + 6);
    C(xC + 6, yC + 7) = A(xA + 6, yA + 0) * B(xB + 0, yB + 7) + A(xA + 6, yA + 1) * B(xB + 1, yB + 7) + A(xA + 6, yA + 2) * B(xB + 2, yB + 7) + A(xA + 6, yA + 3) * B(xB + 3, yB + 7) + A(xA + 6, yA + 4) * B(xB + 4, yB + 7) + A(xA + 6, yA + 5) * B(xB + 5, yB + 7) + A(xA + 6, yA + 6) * B(xB + 6, yB + 7) + A(xA + 6, yA + 7) * B(xB + 7, yB + 7);
    C(xC + 7, yC + 0) = A(xA + 7, yA + 0) * B(xB + 0, yB + 0) + A(xA + 7, yA + 1) * B(xB + 1, yB + 0) + A(xA + 7, yA + 2) * B(xB + 2, yB + 0) + A(xA + 7, yA + 3) * B(xB + 3, yB + 0) + A(xA + 7, yA + 4) * B(xB + 4, yB + 0) + A(xA + 7, yA + 5) * B(xB + 5, yB + 0) + A(xA + 7, yA + 6) * B(xB + 6, yB + 0) + A(xA + 7, yA + 7) * B(xB + 7, yB + 0);
    C(xC + 7, yC + 1) = A(xA + 7, yA + 0) * B(xB + 0, yB + 1) + A(xA + 7, yA + 1) * B(xB + 1, yB + 1) + A(xA + 7, yA + 2) * B(xB + 2, yB + 1) + A(xA + 7, yA + 3) * B(xB + 3, yB + 1) + A(xA + 7, yA + 4) * B(xB + 4, yB + 1) + A(xA + 7, yA + 5) * B(xB + 5, yB + 1) + A(xA + 7, yA + 6) * B(xB + 6, yB + 1) + A(xA + 7, yA + 7) * B(xB + 7, yB + 1);
    C(xC + 7, yC + 2) = A(xA + 7, yA + 0) * B(xB + 0, yB + 2) + A(xA + 7, yA + 1) * B(xB + 1, yB + 2) + A(xA + 7, yA + 2) * B(xB + 2, yB + 2) + A(xA + 7, yA + 3) * B(xB + 3, yB + 2) + A(xA + 7, yA + 4) * B(xB + 4, yB + 2) + A(xA + 7, yA + 5) * B(xB + 5, yB + 2) + A(xA + 7, yA + 6) * B(xB + 6, yB + 2) + A(xA + 7, yA + 7) * B(xB + 7, yB + 2);
    C(xC + 7, yC + 3) = A(xA + 7, yA + 0) * B(xB + 0, yB + 3) + A(xA + 7, yA + 1) * B(xB + 1, yB + 3) + A(xA + 7, yA + 2) * B(xB + 2, yB + 3) + A(xA + 7, yA + 3) * B(xB + 3, yB + 3) + A(xA + 7, yA + 4) * B(xB + 4, yB + 3) + A(xA + 7, yA + 5) * B(xB + 5, yB + 3) + A(xA + 7, yA + 6) * B(xB + 6, yB + 3) + A(xA + 7, yA + 7) * B(xB + 7, yB + 3);
    C(xC + 7, yC + 4) = A(xA + 7, yA + 0) * B(xB + 0, yB + 4) + A(xA + 7, yA + 1) * B(xB + 1, yB + 4) + A(xA + 7, yA + 2) * B(xB + 2, yB + 4) + A(xA + 7, yA + 3) * B(xB + 3, yB + 4) + A(xA + 7, yA + 4) * B(xB + 4, yB + 4) + A(xA + 7, yA + 5) * B(xB + 5, yB + 4) + A(xA + 7, yA + 6) * B(xB + 6, yB + 4) + A(xA + 7, yA + 7) * B(xB + 7, yB + 4);
    C(xC + 7, yC + 5) = A(xA + 7, yA + 0) * B(xB + 0, yB + 5) + A(xA + 7, yA + 1) * B(xB + 1, yB + 5) + A(xA + 7, yA + 2) * B(xB + 2, yB + 5) + A(xA + 7, yA + 3) * B(xB + 3, yB + 5) + A(xA + 7, yA + 4) * B(xB + 4, yB + 5) + A(xA + 7, yA + 5) * B(xB + 5, yB + 5) + A(xA + 7, yA + 6) * B(xB + 6, yB + 5) + A(xA + 7, yA + 7) * B(xB + 7, yB + 5);
    C(xC + 7, yC + 6) = A(xA + 7, yA + 0) * B(xB + 0, yB + 6) + A(xA + 7, yA + 1) * B(xB + 1, yB + 6) + A(xA + 7, yA + 2) * B(xB + 2, yB + 6) + A(xA + 7, yA + 3) * B(xB + 3, yB + 6) + A(xA + 7, yA + 4) * B(xB + 4, yB + 6) + A(xA + 7, yA + 5) * B(xB + 5, yB + 6) + A(xA + 7, yA + 6) * B(xB + 6, yB + 6) + A(xA + 7, yA + 7) * B(xB + 7, yB + 6);
    C(xC + 7, yC + 7) = A(xA + 7, yA + 0) * B(xB + 0, yB + 7) + A(xA + 7, yA + 1) * B(xB + 1, yB + 7) + A(xA + 7, yA + 2) * B(xB + 2, yB + 7) + A(xA + 7, yA + 3) * B(xB + 3, yB + 7) + A(xA + 7, yA + 4) * B(xB + 4, yB + 7) + A(xA + 7, yA + 5) * B(xB + 5, yB + 7) + A(xA + 7, yA + 6) * B(xB + 6, yB + 7) + A(xA + 7, yA + 7) * B(xB + 7, yB + 7);
}

template <class T>
inline void Matrix<T>::base_case4(const Matrix<T>& paddedA, const Matrix<T>& paddedB, Matrix<T>& paddedProd)
{
    paddedProd(0, 0) = paddedA(0, 0) * paddedB(0, 0) + paddedA(0, 1) * paddedB(1, 0) + paddedA(0, 2) * paddedB(2, 0) + paddedA(0, 3) * paddedB(3, 0);
    paddedProd(0, 1) = paddedA(0, 0) * paddedB(0, 1) + paddedA(0, 1) * paddedB(1, 1) + paddedA(0, 2) * paddedB(2, 1) + paddedA(0, 3) * paddedB(3, 1);
    paddedProd(0, 2) = paddedA(0, 0) * paddedB(0, 2) + paddedA(0, 1) * paddedB(1, 2) + paddedA(0, 2) * paddedB(2, 2) + paddedA(0, 3) * paddedB(3, 2);
    paddedProd(0, 3) = paddedA(0, 0) * paddedB(0, 3) + paddedA(0, 1) * paddedB(1, 3) + paddedA(0, 2) * paddedB(2, 3) + paddedA(0, 3) * paddedB(3, 3);
    paddedProd(1, 0) = paddedA(1, 0) * paddedB(0, 0) + paddedA(1, 1) * paddedB(1, 0) + paddedA(1, 2) * paddedB(2, 0) + paddedA(1, 3) * paddedB(3, 0);
    paddedProd(1, 1) = paddedA(1, 0) * paddedB(0, 1) + paddedA(1, 1) * paddedB(1, 1) + paddedA(1, 2) * paddedB(2, 1) + paddedA(1, 3) * paddedB(3, 1);
    paddedProd(1, 2) = paddedA(1, 0) * paddedB(0, 2) + paddedA(1, 1) * paddedB(1, 2) + paddedA(1, 2) * paddedB(2, 2) + paddedA(1, 3) * paddedB(3, 2);
    paddedProd(1, 3) = paddedA(1, 0) * paddedB(0, 3) + paddedA(1, 1) * paddedB(1, 3) + paddedA(1, 2) * paddedB(2, 3) + paddedA(1, 3) * paddedB(3, 3);
    paddedProd(2, 0) = paddedA(2, 0) * paddedB(0, 0) + paddedA(2, 1) * paddedB(1, 0) + paddedA(2, 2) * paddedB(2, 0) + paddedA(2, 3) * paddedB(3, 0);
    paddedProd(2, 1) = paddedA(2, 0) * paddedB(0, 1) + paddedA(2, 1) * paddedB(1, 1) + paddedA(2, 2) * paddedB(2, 1) + paddedA(2, 3) * paddedB(3, 1);
    paddedProd(2, 2) = paddedA(2, 0) * paddedB(0, 2) + paddedA(2, 1) * paddedB(1, 2) + paddedA(2, 2) * paddedB(2, 2) + paddedA(2, 3) * paddedB(3, 2);
    paddedProd(2, 3) = paddedA(2, 0) * paddedB(0, 3) + paddedA(2, 1) * paddedB(1, 3) + paddedA(2, 2) * paddedB(2, 3) + paddedA(2, 3) * paddedB(3, 3);
    paddedProd(3, 0) = paddedA(3, 0) * paddedB(0, 0) + paddedA(3, 1) * paddedB(1, 0) + paddedA(3, 2) * paddedB(2, 0) + paddedA(3, 3) * paddedB(3, 0);
    paddedProd(3, 1) = paddedA(3, 0) * paddedB(0, 1) + paddedA(3, 1) * paddedB(1, 1) + paddedA(3, 2) * paddedB(2, 1) + paddedA(3, 3) * paddedB(3, 1);
    paddedProd(3, 2) = paddedA(3, 0) * paddedB(0, 2) + paddedA(3, 1) * paddedB(1, 2) + paddedA(3, 2) * paddedB(2, 2) + paddedA(3, 3) * paddedB(3, 2);
    paddedProd(3, 3) = paddedA(3, 0) * paddedB(0, 3) + paddedA(3, 1) * paddedB(1, 3) + paddedA(3, 2) * paddedB(2, 3) + paddedA(3, 3) * paddedB(3, 3);
}
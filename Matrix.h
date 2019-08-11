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
    mutable T determinant = 0; // irrelevant in a non-square matrix
    mutable unsigned int rank = 0;
    mutable bool det_updated = true;  // dirty bit, may still be true in the irrelevant case of a non-square matrix
    mutable bool rank_updated = true;  // dirty bit
    std::vector<vector_t> data;       // the matrix's entries, stored as a vector of rows of the matrix, in row major fashion.

    // if matrix is square, the function computes the inverse if exists and updates the rank and det
    // if the matrix is not square, it only updates the rank
    Matrix gaussian_elimination() const
    {
        rank = cols;                // later will perform Gaussian elimination and substract 1 for each column without leading element
        determinant = 1;            // row operations might modify this this value
        Matrix copy(*this);         // the matrix to perform elimination on
        Matrix inverse(rows, rows); // return value
        // set inverse to unit matrix. Later will be eliminated in parallel
        for (unsigned int i = 0; i < rows; ++i) {
            inverse(i, i) = 1;
        }
        unsigned int first_non_zero = 0;                                                                                       // column number in which we expect to find the first nonzero element of some row
        for (unsigned int row = 0; row < rows && first_non_zero < cols;) {                                            // operate on rows one by one
            for (unsigned int lower_row = row + 1; lower_row < rows && copy(row, first_non_zero) == 0; ++lower_row) { // if starts with 0, swap it
                if (copy(lower_row, first_non_zero) != 0) {                                                           // swap the rows and update the data.
                    copy.swap_rows(row, lower_row);
                    inverse.swap_rows(row, lower_row);
                    determinant *= -1;
                }
            }

            if (copy(row, first_non_zero) != 0) {
                T scalar = copy(row, first_non_zero);            // first element
                copy.multiply_row_by_scalar(row, 1 / scalar);    // start with 1.
                inverse.multiply_row_by_scalar(row, 1 / scalar); // parallel elimination.
                determinant *= scalar;
                // set all other rows to zero in the first_non_zero-th index
                for (unsigned int other_row = 0; other_row < rows; ++other_row) {
                    scalar = copy(other_row, first_non_zero);
                    if (other_row != row && scalar != 0) {
                        copy.add_multiplied_row(other_row, row, scalar * (-1));
                        inverse.add_multiplied_row(other_row, row, scalar * (-1));
                    }
                }
                ++first_non_zero; // next
                ++row;            // next
            }
            else {                // no row starting (in first_non_zero) with a nonzero element was found
                determinant = 0;  // not invertible
                ++first_non_zero; // go on
                --rank;
            }
        }

        rank_updated = true;
        det_updated = true;
        if (determinant != 0) {                    // invertible
            inverse.determinant = 1 / determinant; // of course
            inverse.det_updated = true;
            inverse.rank = rows;
            inverse.rank_updated = true;
        }
        return std::move(inverse);
    }

    // max(s: 2^s<=x) x must be nonegative
    unsigned int round_up_power_2(unsigned int x)
    {
        unsigned int i = 0;
        while (x > (1 << i)) {
            ++i;
        }
        return i;
    }

    // set A1 = A1 + B1, where A1 and B1 are square submatrices of A and B, with given starting indices in A, B and the size
    static void add_to_place(const Matrix& A, const Matrix& B, Matrix& C, unsigned int xA, unsigned int yA, unsigned int xB, unsigned int yB, unsigned int xC, unsigned int yC, unsigned int size)
    {
        for (unsigned int i = 0; i < size; i++) {
            for (unsigned int j = 0; j < size; j++) {
                C(xC + i, yC + j) = A(xA + i, yA + j) + B(xB + i, yB + j);
            }
        }
    }

    // set A1 = A1 - B1, where A1 and B1 are square submatrices of A and B, with given starting indices in A, B and the size
    static void sub_to_place(const Matrix& A, const Matrix& B, Matrix& C, unsigned int xA, unsigned int yA, unsigned int xB, unsigned int yB, unsigned int xC, unsigned int yC, unsigned int size)
    {
        for (unsigned int i = 0; i < size; i++) {
            for (unsigned int j = 0; j < size; j++) {
                C(xC + i, yC + j) = A(xA + i, yA + j) - B(xB + i, yB + j);
            }
        }
    }

    // set A1 = B1, where A1 and B1 are square submatrices of A and B, with given starting indices in A, B and the size
    static void copy(Matrix& A, const Matrix& B, unsigned int xA, unsigned int yA, unsigned int xB, unsigned int yB, unsigned int size)
    {
        for (unsigned int i = 0; i < size; i++) {
            for (unsigned int j = 0; j < size; j++) {
                A(xA + i, yA + j) = B(xB + i, yB + j);
            }
        }
    }

    // function just for Strassen multiplication
    // loop unrolling code for multiplying two 8X8 matrices
    static void base_case8(const Matrix& A, const Matrix& B, Matrix& C, unsigned int xA, unsigned int yA, unsigned int xB, unsigned int yB, unsigned int xC, unsigned int yC);

    // function just for Strassen multiplication
    // loop unrolling code for multiplying two 4X4 matrices
    static void base_case4(const Matrix& paddedA, const Matrix& paddedB, Matrix& paddedProd);

    // function just for Strassen multiplication
    // recursive multiplication the of square matrices of dimension = power of two >= 8
    // matrices may be given as submatrices of larger matrices A and B, by starting indices and actual size
    static void add_mul(const Matrix& A, const Matrix& B, Matrix& C, unsigned int xA, unsigned int yA, unsigned int xB, unsigned int yB, unsigned int xC, unsigned int yC, unsigned int size)
    {
        // the actual Strassen login: partition to 4 blocks of half size, and recursive operation:

        if (size == 8) { // base case
            Matrix::base_case8(A, B, C, xA, yA, xB, yB, xC, yC);
            return;
        }

        unsigned int half_size = size / 2; // blocks size

        // blocks which are meant to store all the relevant information for the calculations
        // of the new blocks that will compose the result.
        Matrix X(half_size, half_size);
        Matrix Y(half_size, half_size);

        // in Strassen's algorithm denote:
        // M1 = (A11 + A22) * (B11 + B22)
        // M2 = (A21 + A22) * B11
        // M3 = A11 * (B12 - B22)
        // M4 = A22 * (B21 - B11)
        // M5 = (A11 + A12) * B22
        // M6 = (A21 - A11) * (B11 + B12)
        // M7 = (A12 - A22) * (B21 + B22)
        //
        // and then:
        //
        // C11 = M1 + M4 - M5 + M7
        // C12 = M3 + M5
        // C21 = M2 + M4
        // C22 = M1 - M2 + M3 + M6
        //
        // where
        //
        //     A11 | A12S
        // A = ----------
        //     A21 | A22
        //
        // and similarly to B and C.
        // here we don't initialize Aij, Bij and Mi, but instead perform most of them in the
        // designated cells: X, Y and M.

        // M1
        add_to_place(A, A, X, xA, yA, xA + half_size, yA + half_size, 0, 0, half_size);
        add_to_place(B, B, Y, xB, yB, xB + half_size, yB + half_size, 0, 0, half_size);

        sub_to_place(C, C, C, xC + half_size, yC + half_size, xC, yC, xC + half_size, yC + half_size, half_size);
        add_mul(X, Y, C, 0, 0, 0, 0, xC, yC, half_size);
        add_to_place(C, C, C, xC + half_size, yC + half_size, xC, yC, xC + half_size, yC + half_size, half_size);

        // M2
        add_to_place(A, A, X, xA + half_size, yA, xA + half_size, yA + half_size, 0, 0, half_size);

        Y *= 0;
        add_mul(X, B, Y, 0, 0, xB, yB, 0, 0, half_size);
        add_to_place(C, Y, C, xC + half_size, yC, 0, 0, xC + half_size, yC, half_size);
        sub_to_place(C, Y, C, xC + half_size, yC + half_size, 0, 0, xC + half_size, yC + half_size, half_size);

        // M3
        sub_to_place(B, B, Y, xB, yB + half_size, xB + half_size, yB + half_size, 0, 0, half_size);

        X *= 0;
        add_mul(A, Y, X, xA, yA, 0, 0, 0, 0, half_size);
        add_to_place(C, X, C, xC, yC + half_size, 0, 0, xC, yC + half_size, half_size);
        add_to_place(C, X, C, xC + half_size, yC + half_size, 0, 0, xC + half_size, yC + half_size, half_size);

        // M4
        sub_to_place(B, B, Y, xB + half_size, yB, xB, yB, 0, 0, half_size);

        X *= 0;
        add_mul(A, Y, X, xA + half_size, yA + half_size, 0, 0, 0, 0, half_size);
        add_to_place(C, X, C, xC, yC, 0, 0, xC, yC, half_size);
        add_to_place(C, X, C, xC + half_size, yC, 0, 0, xC + half_size, yC, half_size);

        // M5
        add_to_place(A, A, X, xA, yA, xA, yA + half_size, 0, 0, half_size);

        Y *= 0;
        add_mul(X, B, Y, 0, 0, xB + half_size, yB + half_size, 0, 0, half_size);
        sub_to_place(C, Y, C, xC, yC, 0, 0, xC, yC, half_size);
        add_to_place(C, Y, C, xC, yC + half_size, 0, 0, xC, yC + half_size, half_size);

        // M6
        sub_to_place(A, A, X, xA + half_size, yA, xA, yA, 0, 0, half_size);
        add_to_place(B, B, Y, xB, yB, xB, yB + half_size, 0, 0, half_size);

        add_mul(X, Y, C, 0, 0, 0, 0, xC + half_size, yC + half_size, half_size);

        // M7
        sub_to_place(A, A, X, xA, yA + half_size, xA + half_size, yA + half_size, 0, 0, half_size);
        add_to_place(B, B, Y, xB + half_size, yB, xB + half_size, yB + half_size, 0, 0, half_size);

        add_mul(X, Y, C, 0, 0, 0, 0, xC, yC, half_size);
    }

public:
    // default constructor. Constructs a 1X1 with 0
    explicit Matrix() : rows(1), cols(1), data(std::vector<vector_t>(1, vector_t(1, 0)))
    {
    }

    // constructs an rXc matrix, filled with 0's
    explicit Matrix(unsigned int r, unsigned int c) : rows(r), cols(c), data(std::vector<vector_t>(rows, vector_t(cols, 0)))
    {
    }

    // constructor by const reference to std::vector of vector_t (the entries)
    // assumes the std::vector's in entries are of the same length
    explicit Matrix(const std::vector<vector_t>& entries) : rows(entries.size()), cols(entries[0].size()), data(entries), det_updated(false), rank_updated(false)
    {
    }

    // constructor by rvalue reference to std::vector of vector_t (the entries)
    // assumes the std::vector's in entries are of the same length
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

    // move assignment by rvalue reference
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

    // return whether the given matrix has the same entries.
    // for programmers: it also updates the det or rank if the matrices are equal and one of the has the data up to date.
    bool operator==(const Matrix& m) const
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

    // return whether the given matrix has the different entries.
    bool operator!=(const Matrix& m) const
    {
        return !(*this == m);
    }

    // returns the number of rows
    unsigned int get_rows() const
    {
        return rows;
    }

    // returns the number of columns
    unsigned int get_cols() const
    {
        return cols;
    }

    // get: matrix(i, j) - the matrix's entry at the i'th row and j'th column, zero based
    // indices must be in the appropriate range.
    const T& operator()(unsigned int i, unsigned int j) const
    {
        return data[i][j];
    }

    // set: matrix(i, j) - the matrix's entry at the i'th row and j'th column, zero based
    // indices must be in the appropriate range.
    T& operator()(unsigned int i, unsigned int j)
    {
        rank_updated = false; // entries might change now
        det_updated = false;
        return data[i][j];
    }

    // get: matrix(i) is the matrix's i'th row, zero based
    // index should be in the appropriate range.
    const vector_t& operator()(unsigned int i) const
    {
        return data[i];
    }

    // set: matrix(i) is the matrix's i'th row, zero based
    // index should be in the appropriate range.
    vector_t& operator()(unsigned int i)
    {
        det_updated = false; // entries might change now
        rank_updated = false;
        return data[i];
    }

    // changes the matrix to its transposed form
    // assumes the matrix is squared!
    void transpose_in_place()
    {
        for (unsigned int i = 0; i < rows; ++i) {
            for (unsigned int j = 0; j < i; ++j) {
                std::swap((*this)(i, j), (*this)(j, i));
            }
        }
    }

    // returns the transposed matrix.
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
    // assumes they are the same size.
    Matrix& operator+=(const Matrix& m)
    {
        for (unsigned int i = 0; i < rows; ++i) {
            for (unsigned int j = 0; j < cols; ++j) {
                (*this)(i, j) += m(i, j);
            }
        }
        rank_updated = false; // will take time to update - might hurt performance when irrelevant
        det_updated = false;
        return *this;
    }

    // matrices addition
    // assumes they are the same size.
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
    // assumes they are mutiplicable.
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

    // faster multiplication of two matrices.
    // assumes they are mutiplicable.
    Matrix strassen(const Matrix<T>& B)
    {
        Matrix A = *this;
        unsigned int power = round_up_power_2(std::max({A.get_rows(), A.get_cols(), B.get_cols()}));
        unsigned int size = 1 << power;

        // padding the matrices by zeroes to be parts of a power of 2 dimensional squared matrices.
        Matrix paddedA(size, size);
        Matrix paddedB(size, size);

        for (unsigned int i = 0; i < A.get_rows(); ++i) {
            for (unsigned int j = 0; j < A.get_cols(); ++j) {
                paddedA(i, j) = A(i, j);
            }
        }

        for (unsigned int i = 0; i < B.get_rows(); ++i) {
            for (unsigned int j = 0; j < B.get_cols(); ++j) {
                paddedB(i, j) = B(i, j);
            }
        }

        Matrix paddedProd = Matrix(size, size);
        if (size == 1) { // case too small to mul
            paddedProd(0, 0) = A(0, 0) * B(0, 0);
        }

        else if (size == 2) { // case too small to mul
            paddedProd(0, 0) = A(0, 0) * B(0, 0) + A(0, 1) * B(1, 0);
            paddedProd(0, 1) = A(0, 0) * B(0, 1) + A(0, 1) * B(1, 1);
            paddedProd(1, 0) = A(1, 0) * B(0, 0) + A(1, 1) * B(1, 0);
            paddedProd(1, 1) = A(1, 0) * B(0, 1) + A(1, 1) * B(1, 1);
        }

        else if (size == 4) { // case too small to mul
            base_case4(paddedA, paddedB, paddedProd);
        }

        else {
            // calculate the product with add_mul()
            add_mul(paddedA, paddedB, paddedProd, 0, 0, 0, 0, 0, 0, size);
        }

        // unpad the product
        Matrix prod = Matrix(A.get_rows(), B.get_cols());
        for (unsigned int i = 0; i < prod.get_rows(); ++i) {
            for (unsigned int j = 0; j < prod.get_cols(); ++j) {
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
    // assumes they are the same size.
    Matrix operator-() const
    {
        return (*this) * (-1);
    }

    // in-place matrices substraction
    // assumes they are the same size.
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
    // assumes they are the same size.
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
    // assumes they are multiplicable.
    Matrix& operator*=(const T& s)
    {
        bool det_up = det_updated;   // assignment might mark as dirty
        bool rank_up = rank_updated; // assignment might mark as dirty

        for (unsigned int i = 0; i < rows; ++i) {
            for (unsigned int j = 0; j < cols; ++j) {
                (*this)(i, j) *= s;
            }
        }
        det_updated = det_up;
        rank_updated = rank_up;

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
    // assumes they are multiplicable.
    Vector<T> operator*(const Vector<T>& v) const
    {
        Vector<T> res(this->get_rows()); // output
        for (unsigned int i = 0; i < this->get_rows(); ++i) {
            res(i) = 0;
            for (unsigned int j = 0; j < this->get_cols(); ++j) {
                res(i) += (*this)(i, j) * v(j);
            }
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

    // gets two rows and swap them
    // the rows should have a valid range
    void swap_rows(unsigned int i, unsigned int j)
    {
        // remember the values so that the "set" operations will not mark them as dirty
        bool det_up = det_updated;
        bool rank_up = rank_updated;

        std::swap((*this)(i), (*this)(j));
        determinant *= -1;

        det_updated = det_up;
        rank_updated = rank_up;
    }

    // gets a row of the matrix and a scalar, and multiplies (in place) the row by the scalar.
    void multiply_row_by_scalar(unsigned int row, T scalar)
    {
        // remember the values so that the "set" operations will not mark them as dirty
        bool det_up = det_updated;
        bool rank_up = rank_updated;

        for (unsigned int i = 0; i < cols; ++i) { // multiply the row
            (*this)(row, i) *= scalar;
        }
        determinant *= scalar;

        det_updated = det_up;
        if (scalar == 0) {
            rank_updated = false;
        }
        else {
            rank_updated = rank_up;
        }
    }

    // adds a multiplication of some row to another row (in place)
    void add_multiplied_row(unsigned int row1, unsigned int row2, T scalar)
    {
        // remember the values so that the "set" operations will not mark them as dirty
        bool det_up = det_updated;
        bool rank_up = rank_updated;

        for (unsigned int i = 0; i < cols; ++i) {
            (*this)(row1, i) += (*this)(row2, i) * scalar;
        }

        det_updated = det_up;
        rank_updated = rank_up;
    }

    // returns the determinant of a matrix
    // assumes the matrix is square
    T get_det() const
    {
        if (!det_updated) {
            gaussian_elimination();
        }
        return determinant;
    }
    // assumes the matrix has inverse
    Matrix find_inverse() const
    {
        return gaussian_elimination();
    }
    //
    unsigned int get_rank() const
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
    strm << "["; // formatting nicely
    for (unsigned int i = 0; i < m.get_rows(); ++i) {
        strm << "("; // row opener
        for (unsigned int j = 0; j < m.get_cols() - 1; ++j) {
            strm << m(i, j) << ", "; // values seperator
        }
        strm << m(i, m.get_cols() - 1) << ")\n";
    }
    strm << "]" << std::endl; // end
    return strm;
}

template <class T>
inline void Matrix<T>::base_case8(const Matrix<T>& A, const Matrix<T>& B, Matrix<T>& C, unsigned int xA, unsigned int yA, unsigned int xB, unsigned int yB, unsigned int xC, unsigned int yC)
{
    C(xC + 0, yC + 0) += A(xA + 0, yA + 0) * B(xB + 0, yB + 0) + A(xA + 0, yA + 1) * B(xB + 1, yB + 0) + A(xA + 0, yA + 2) * B(xB + 2, yB + 0) + A(xA + 0, yA + 3) * B(xB + 3, yB + 0) + A(xA + 0, yA + 4) * B(xB + 4, yB + 0) + A(xA + 0, yA + 5) * B(xB + 5, yB + 0) + A(xA + 0, yA + 6) * B(xB + 6, yB + 0) + A(xA + 0, yA + 7) * B(xB + 7, yB + 0);
    C(xC + 0, yC + 1) += A(xA + 0, yA + 0) * B(xB + 0, yB + 1) + A(xA + 0, yA + 1) * B(xB + 1, yB + 1) + A(xA + 0, yA + 2) * B(xB + 2, yB + 1) + A(xA + 0, yA + 3) * B(xB + 3, yB + 1) + A(xA + 0, yA + 4) * B(xB + 4, yB + 1) + A(xA + 0, yA + 5) * B(xB + 5, yB + 1) + A(xA + 0, yA + 6) * B(xB + 6, yB + 1) + A(xA + 0, yA + 7) * B(xB + 7, yB + 1);
    C(xC + 0, yC + 2) += A(xA + 0, yA + 0) * B(xB + 0, yB + 2) + A(xA + 0, yA + 1) * B(xB + 1, yB + 2) + A(xA + 0, yA + 2) * B(xB + 2, yB + 2) + A(xA + 0, yA + 3) * B(xB + 3, yB + 2) + A(xA + 0, yA + 4) * B(xB + 4, yB + 2) + A(xA + 0, yA + 5) * B(xB + 5, yB + 2) + A(xA + 0, yA + 6) * B(xB + 6, yB + 2) + A(xA + 0, yA + 7) * B(xB + 7, yB + 2);
    C(xC + 0, yC + 3) += A(xA + 0, yA + 0) * B(xB + 0, yB + 3) + A(xA + 0, yA + 1) * B(xB + 1, yB + 3) + A(xA + 0, yA + 2) * B(xB + 2, yB + 3) + A(xA + 0, yA + 3) * B(xB + 3, yB + 3) + A(xA + 0, yA + 4) * B(xB + 4, yB + 3) + A(xA + 0, yA + 5) * B(xB + 5, yB + 3) + A(xA + 0, yA + 6) * B(xB + 6, yB + 3) + A(xA + 0, yA + 7) * B(xB + 7, yB + 3);
    C(xC + 0, yC + 4) += A(xA + 0, yA + 0) * B(xB + 0, yB + 4) + A(xA + 0, yA + 1) * B(xB + 1, yB + 4) + A(xA + 0, yA + 2) * B(xB + 2, yB + 4) + A(xA + 0, yA + 3) * B(xB + 3, yB + 4) + A(xA + 0, yA + 4) * B(xB + 4, yB + 4) + A(xA + 0, yA + 5) * B(xB + 5, yB + 4) + A(xA + 0, yA + 6) * B(xB + 6, yB + 4) + A(xA + 0, yA + 7) * B(xB + 7, yB + 4);
    C(xC + 0, yC + 5) += A(xA + 0, yA + 0) * B(xB + 0, yB + 5) + A(xA + 0, yA + 1) * B(xB + 1, yB + 5) + A(xA + 0, yA + 2) * B(xB + 2, yB + 5) + A(xA + 0, yA + 3) * B(xB + 3, yB + 5) + A(xA + 0, yA + 4) * B(xB + 4, yB + 5) + A(xA + 0, yA + 5) * B(xB + 5, yB + 5) + A(xA + 0, yA + 6) * B(xB + 6, yB + 5) + A(xA + 0, yA + 7) * B(xB + 7, yB + 5);
    C(xC + 0, yC + 6) += A(xA + 0, yA + 0) * B(xB + 0, yB + 6) + A(xA + 0, yA + 1) * B(xB + 1, yB + 6) + A(xA + 0, yA + 2) * B(xB + 2, yB + 6) + A(xA + 0, yA + 3) * B(xB + 3, yB + 6) + A(xA + 0, yA + 4) * B(xB + 4, yB + 6) + A(xA + 0, yA + 5) * B(xB + 5, yB + 6) + A(xA + 0, yA + 6) * B(xB + 6, yB + 6) + A(xA + 0, yA + 7) * B(xB + 7, yB + 6);
    C(xC + 0, yC + 7) += A(xA + 0, yA + 0) * B(xB + 0, yB + 7) + A(xA + 0, yA + 1) * B(xB + 1, yB + 7) + A(xA + 0, yA + 2) * B(xB + 2, yB + 7) + A(xA + 0, yA + 3) * B(xB + 3, yB + 7) + A(xA + 0, yA + 4) * B(xB + 4, yB + 7) + A(xA + 0, yA + 5) * B(xB + 5, yB + 7) + A(xA + 0, yA + 6) * B(xB + 6, yB + 7) + A(xA + 0, yA + 7) * B(xB + 7, yB + 7);
    C(xC + 1, yC + 0) += A(xA + 1, yA + 0) * B(xB + 0, yB + 0) + A(xA + 1, yA + 1) * B(xB + 1, yB + 0) + A(xA + 1, yA + 2) * B(xB + 2, yB + 0) + A(xA + 1, yA + 3) * B(xB + 3, yB + 0) + A(xA + 1, yA + 4) * B(xB + 4, yB + 0) + A(xA + 1, yA + 5) * B(xB + 5, yB + 0) + A(xA + 1, yA + 6) * B(xB + 6, yB + 0) + A(xA + 1, yA + 7) * B(xB + 7, yB + 0);
    C(xC + 1, yC + 1) += A(xA + 1, yA + 0) * B(xB + 0, yB + 1) + A(xA + 1, yA + 1) * B(xB + 1, yB + 1) + A(xA + 1, yA + 2) * B(xB + 2, yB + 1) + A(xA + 1, yA + 3) * B(xB + 3, yB + 1) + A(xA + 1, yA + 4) * B(xB + 4, yB + 1) + A(xA + 1, yA + 5) * B(xB + 5, yB + 1) + A(xA + 1, yA + 6) * B(xB + 6, yB + 1) + A(xA + 1, yA + 7) * B(xB + 7, yB + 1);
    C(xC + 1, yC + 2) += A(xA + 1, yA + 0) * B(xB + 0, yB + 2) + A(xA + 1, yA + 1) * B(xB + 1, yB + 2) + A(xA + 1, yA + 2) * B(xB + 2, yB + 2) + A(xA + 1, yA + 3) * B(xB + 3, yB + 2) + A(xA + 1, yA + 4) * B(xB + 4, yB + 2) + A(xA + 1, yA + 5) * B(xB + 5, yB + 2) + A(xA + 1, yA + 6) * B(xB + 6, yB + 2) + A(xA + 1, yA + 7) * B(xB + 7, yB + 2);
    C(xC + 1, yC + 3) += A(xA + 1, yA + 0) * B(xB + 0, yB + 3) + A(xA + 1, yA + 1) * B(xB + 1, yB + 3) + A(xA + 1, yA + 2) * B(xB + 2, yB + 3) + A(xA + 1, yA + 3) * B(xB + 3, yB + 3) + A(xA + 1, yA + 4) * B(xB + 4, yB + 3) + A(xA + 1, yA + 5) * B(xB + 5, yB + 3) + A(xA + 1, yA + 6) * B(xB + 6, yB + 3) + A(xA + 1, yA + 7) * B(xB + 7, yB + 3);
    C(xC + 1, yC + 4) += A(xA + 1, yA + 0) * B(xB + 0, yB + 4) + A(xA + 1, yA + 1) * B(xB + 1, yB + 4) + A(xA + 1, yA + 2) * B(xB + 2, yB + 4) + A(xA + 1, yA + 3) * B(xB + 3, yB + 4) + A(xA + 1, yA + 4) * B(xB + 4, yB + 4) + A(xA + 1, yA + 5) * B(xB + 5, yB + 4) + A(xA + 1, yA + 6) * B(xB + 6, yB + 4) + A(xA + 1, yA + 7) * B(xB + 7, yB + 4);
    C(xC + 1, yC + 5) += A(xA + 1, yA + 0) * B(xB + 0, yB + 5) + A(xA + 1, yA + 1) * B(xB + 1, yB + 5) + A(xA + 1, yA + 2) * B(xB + 2, yB + 5) + A(xA + 1, yA + 3) * B(xB + 3, yB + 5) + A(xA + 1, yA + 4) * B(xB + 4, yB + 5) + A(xA + 1, yA + 5) * B(xB + 5, yB + 5) + A(xA + 1, yA + 6) * B(xB + 6, yB + 5) + A(xA + 1, yA + 7) * B(xB + 7, yB + 5);
    C(xC + 1, yC + 6) += A(xA + 1, yA + 0) * B(xB + 0, yB + 6) + A(xA + 1, yA + 1) * B(xB + 1, yB + 6) + A(xA + 1, yA + 2) * B(xB + 2, yB + 6) + A(xA + 1, yA + 3) * B(xB + 3, yB + 6) + A(xA + 1, yA + 4) * B(xB + 4, yB + 6) + A(xA + 1, yA + 5) * B(xB + 5, yB + 6) + A(xA + 1, yA + 6) * B(xB + 6, yB + 6) + A(xA + 1, yA + 7) * B(xB + 7, yB + 6);
    C(xC + 1, yC + 7) += A(xA + 1, yA + 0) * B(xB + 0, yB + 7) + A(xA + 1, yA + 1) * B(xB + 1, yB + 7) + A(xA + 1, yA + 2) * B(xB + 2, yB + 7) + A(xA + 1, yA + 3) * B(xB + 3, yB + 7) + A(xA + 1, yA + 4) * B(xB + 4, yB + 7) + A(xA + 1, yA + 5) * B(xB + 5, yB + 7) + A(xA + 1, yA + 6) * B(xB + 6, yB + 7) + A(xA + 1, yA + 7) * B(xB + 7, yB + 7);
    C(xC + 2, yC + 0) += A(xA + 2, yA + 0) * B(xB + 0, yB + 0) + A(xA + 2, yA + 1) * B(xB + 1, yB + 0) + A(xA + 2, yA + 2) * B(xB + 2, yB + 0) + A(xA + 2, yA + 3) * B(xB + 3, yB + 0) + A(xA + 2, yA + 4) * B(xB + 4, yB + 0) + A(xA + 2, yA + 5) * B(xB + 5, yB + 0) + A(xA + 2, yA + 6) * B(xB + 6, yB + 0) + A(xA + 2, yA + 7) * B(xB + 7, yB + 0);
    C(xC + 2, yC + 1) += A(xA + 2, yA + 0) * B(xB + 0, yB + 1) + A(xA + 2, yA + 1) * B(xB + 1, yB + 1) + A(xA + 2, yA + 2) * B(xB + 2, yB + 1) + A(xA + 2, yA + 3) * B(xB + 3, yB + 1) + A(xA + 2, yA + 4) * B(xB + 4, yB + 1) + A(xA + 2, yA + 5) * B(xB + 5, yB + 1) + A(xA + 2, yA + 6) * B(xB + 6, yB + 1) + A(xA + 2, yA + 7) * B(xB + 7, yB + 1);
    C(xC + 2, yC + 2) += A(xA + 2, yA + 0) * B(xB + 0, yB + 2) + A(xA + 2, yA + 1) * B(xB + 1, yB + 2) + A(xA + 2, yA + 2) * B(xB + 2, yB + 2) + A(xA + 2, yA + 3) * B(xB + 3, yB + 2) + A(xA + 2, yA + 4) * B(xB + 4, yB + 2) + A(xA + 2, yA + 5) * B(xB + 5, yB + 2) + A(xA + 2, yA + 6) * B(xB + 6, yB + 2) + A(xA + 2, yA + 7) * B(xB + 7, yB + 2);
    C(xC + 2, yC + 3) += A(xA + 2, yA + 0) * B(xB + 0, yB + 3) + A(xA + 2, yA + 1) * B(xB + 1, yB + 3) + A(xA + 2, yA + 2) * B(xB + 2, yB + 3) + A(xA + 2, yA + 3) * B(xB + 3, yB + 3) + A(xA + 2, yA + 4) * B(xB + 4, yB + 3) + A(xA + 2, yA + 5) * B(xB + 5, yB + 3) + A(xA + 2, yA + 6) * B(xB + 6, yB + 3) + A(xA + 2, yA + 7) * B(xB + 7, yB + 3);
    C(xC + 2, yC + 4) += A(xA + 2, yA + 0) * B(xB + 0, yB + 4) + A(xA + 2, yA + 1) * B(xB + 1, yB + 4) + A(xA + 2, yA + 2) * B(xB + 2, yB + 4) + A(xA + 2, yA + 3) * B(xB + 3, yB + 4) + A(xA + 2, yA + 4) * B(xB + 4, yB + 4) + A(xA + 2, yA + 5) * B(xB + 5, yB + 4) + A(xA + 2, yA + 6) * B(xB + 6, yB + 4) + A(xA + 2, yA + 7) * B(xB + 7, yB + 4);
    C(xC + 2, yC + 5) += A(xA + 2, yA + 0) * B(xB + 0, yB + 5) + A(xA + 2, yA + 1) * B(xB + 1, yB + 5) + A(xA + 2, yA + 2) * B(xB + 2, yB + 5) + A(xA + 2, yA + 3) * B(xB + 3, yB + 5) + A(xA + 2, yA + 4) * B(xB + 4, yB + 5) + A(xA + 2, yA + 5) * B(xB + 5, yB + 5) + A(xA + 2, yA + 6) * B(xB + 6, yB + 5) + A(xA + 2, yA + 7) * B(xB + 7, yB + 5);
    C(xC + 2, yC + 6) += A(xA + 2, yA + 0) * B(xB + 0, yB + 6) + A(xA + 2, yA + 1) * B(xB + 1, yB + 6) + A(xA + 2, yA + 2) * B(xB + 2, yB + 6) + A(xA + 2, yA + 3) * B(xB + 3, yB + 6) + A(xA + 2, yA + 4) * B(xB + 4, yB + 6) + A(xA + 2, yA + 5) * B(xB + 5, yB + 6) + A(xA + 2, yA + 6) * B(xB + 6, yB + 6) + A(xA + 2, yA + 7) * B(xB + 7, yB + 6);
    C(xC + 2, yC + 7) += A(xA + 2, yA + 0) * B(xB + 0, yB + 7) + A(xA + 2, yA + 1) * B(xB + 1, yB + 7) + A(xA + 2, yA + 2) * B(xB + 2, yB + 7) + A(xA + 2, yA + 3) * B(xB + 3, yB + 7) + A(xA + 2, yA + 4) * B(xB + 4, yB + 7) + A(xA + 2, yA + 5) * B(xB + 5, yB + 7) + A(xA + 2, yA + 6) * B(xB + 6, yB + 7) + A(xA + 2, yA + 7) * B(xB + 7, yB + 7);
    C(xC + 3, yC + 0) += A(xA + 3, yA + 0) * B(xB + 0, yB + 0) + A(xA + 3, yA + 1) * B(xB + 1, yB + 0) + A(xA + 3, yA + 2) * B(xB + 2, yB + 0) + A(xA + 3, yA + 3) * B(xB + 3, yB + 0) + A(xA + 3, yA + 4) * B(xB + 4, yB + 0) + A(xA + 3, yA + 5) * B(xB + 5, yB + 0) + A(xA + 3, yA + 6) * B(xB + 6, yB + 0) + A(xA + 3, yA + 7) * B(xB + 7, yB + 0);
    C(xC + 3, yC + 1) += A(xA + 3, yA + 0) * B(xB + 0, yB + 1) + A(xA + 3, yA + 1) * B(xB + 1, yB + 1) + A(xA + 3, yA + 2) * B(xB + 2, yB + 1) + A(xA + 3, yA + 3) * B(xB + 3, yB + 1) + A(xA + 3, yA + 4) * B(xB + 4, yB + 1) + A(xA + 3, yA + 5) * B(xB + 5, yB + 1) + A(xA + 3, yA + 6) * B(xB + 6, yB + 1) + A(xA + 3, yA + 7) * B(xB + 7, yB + 1);
    C(xC + 3, yC + 2) += A(xA + 3, yA + 0) * B(xB + 0, yB + 2) + A(xA + 3, yA + 1) * B(xB + 1, yB + 2) + A(xA + 3, yA + 2) * B(xB + 2, yB + 2) + A(xA + 3, yA + 3) * B(xB + 3, yB + 2) + A(xA + 3, yA + 4) * B(xB + 4, yB + 2) + A(xA + 3, yA + 5) * B(xB + 5, yB + 2) + A(xA + 3, yA + 6) * B(xB + 6, yB + 2) + A(xA + 3, yA + 7) * B(xB + 7, yB + 2);
    C(xC + 3, yC + 3) += A(xA + 3, yA + 0) * B(xB + 0, yB + 3) + A(xA + 3, yA + 1) * B(xB + 1, yB + 3) + A(xA + 3, yA + 2) * B(xB + 2, yB + 3) + A(xA + 3, yA + 3) * B(xB + 3, yB + 3) + A(xA + 3, yA + 4) * B(xB + 4, yB + 3) + A(xA + 3, yA + 5) * B(xB + 5, yB + 3) + A(xA + 3, yA + 6) * B(xB + 6, yB + 3) + A(xA + 3, yA + 7) * B(xB + 7, yB + 3);
    C(xC + 3, yC + 4) += A(xA + 3, yA + 0) * B(xB + 0, yB + 4) + A(xA + 3, yA + 1) * B(xB + 1, yB + 4) + A(xA + 3, yA + 2) * B(xB + 2, yB + 4) + A(xA + 3, yA + 3) * B(xB + 3, yB + 4) + A(xA + 3, yA + 4) * B(xB + 4, yB + 4) + A(xA + 3, yA + 5) * B(xB + 5, yB + 4) + A(xA + 3, yA + 6) * B(xB + 6, yB + 4) + A(xA + 3, yA + 7) * B(xB + 7, yB + 4);
    C(xC + 3, yC + 5) += A(xA + 3, yA + 0) * B(xB + 0, yB + 5) + A(xA + 3, yA + 1) * B(xB + 1, yB + 5) + A(xA + 3, yA + 2) * B(xB + 2, yB + 5) + A(xA + 3, yA + 3) * B(xB + 3, yB + 5) + A(xA + 3, yA + 4) * B(xB + 4, yB + 5) + A(xA + 3, yA + 5) * B(xB + 5, yB + 5) + A(xA + 3, yA + 6) * B(xB + 6, yB + 5) + A(xA + 3, yA + 7) * B(xB + 7, yB + 5);
    C(xC + 3, yC + 6) += A(xA + 3, yA + 0) * B(xB + 0, yB + 6) + A(xA + 3, yA + 1) * B(xB + 1, yB + 6) + A(xA + 3, yA + 2) * B(xB + 2, yB + 6) + A(xA + 3, yA + 3) * B(xB + 3, yB + 6) + A(xA + 3, yA + 4) * B(xB + 4, yB + 6) + A(xA + 3, yA + 5) * B(xB + 5, yB + 6) + A(xA + 3, yA + 6) * B(xB + 6, yB + 6) + A(xA + 3, yA + 7) * B(xB + 7, yB + 6);
    C(xC + 3, yC + 7) += A(xA + 3, yA + 0) * B(xB + 0, yB + 7) + A(xA + 3, yA + 1) * B(xB + 1, yB + 7) + A(xA + 3, yA + 2) * B(xB + 2, yB + 7) + A(xA + 3, yA + 3) * B(xB + 3, yB + 7) + A(xA + 3, yA + 4) * B(xB + 4, yB + 7) + A(xA + 3, yA + 5) * B(xB + 5, yB + 7) + A(xA + 3, yA + 6) * B(xB + 6, yB + 7) + A(xA + 3, yA + 7) * B(xB + 7, yB + 7);
    C(xC + 4, yC + 0) += A(xA + 4, yA + 0) * B(xB + 0, yB + 0) + A(xA + 4, yA + 1) * B(xB + 1, yB + 0) + A(xA + 4, yA + 2) * B(xB + 2, yB + 0) + A(xA + 4, yA + 3) * B(xB + 3, yB + 0) + A(xA + 4, yA + 4) * B(xB + 4, yB + 0) + A(xA + 4, yA + 5) * B(xB + 5, yB + 0) + A(xA + 4, yA + 6) * B(xB + 6, yB + 0) + A(xA + 4, yA + 7) * B(xB + 7, yB + 0);
    C(xC + 4, yC + 1) += A(xA + 4, yA + 0) * B(xB + 0, yB + 1) + A(xA + 4, yA + 1) * B(xB + 1, yB + 1) + A(xA + 4, yA + 2) * B(xB + 2, yB + 1) + A(xA + 4, yA + 3) * B(xB + 3, yB + 1) + A(xA + 4, yA + 4) * B(xB + 4, yB + 1) + A(xA + 4, yA + 5) * B(xB + 5, yB + 1) + A(xA + 4, yA + 6) * B(xB + 6, yB + 1) + A(xA + 4, yA + 7) * B(xB + 7, yB + 1);
    C(xC + 4, yC + 2) += A(xA + 4, yA + 0) * B(xB + 0, yB + 2) + A(xA + 4, yA + 1) * B(xB + 1, yB + 2) + A(xA + 4, yA + 2) * B(xB + 2, yB + 2) + A(xA + 4, yA + 3) * B(xB + 3, yB + 2) + A(xA + 4, yA + 4) * B(xB + 4, yB + 2) + A(xA + 4, yA + 5) * B(xB + 5, yB + 2) + A(xA + 4, yA + 6) * B(xB + 6, yB + 2) + A(xA + 4, yA + 7) * B(xB + 7, yB + 2);
    C(xC + 4, yC + 3) += A(xA + 4, yA + 0) * B(xB + 0, yB + 3) + A(xA + 4, yA + 1) * B(xB + 1, yB + 3) + A(xA + 4, yA + 2) * B(xB + 2, yB + 3) + A(xA + 4, yA + 3) * B(xB + 3, yB + 3) + A(xA + 4, yA + 4) * B(xB + 4, yB + 3) + A(xA + 4, yA + 5) * B(xB + 5, yB + 3) + A(xA + 4, yA + 6) * B(xB + 6, yB + 3) + A(xA + 4, yA + 7) * B(xB + 7, yB + 3);
    C(xC + 4, yC + 4) += A(xA + 4, yA + 0) * B(xB + 0, yB + 4) + A(xA + 4, yA + 1) * B(xB + 1, yB + 4) + A(xA + 4, yA + 2) * B(xB + 2, yB + 4) + A(xA + 4, yA + 3) * B(xB + 3, yB + 4) + A(xA + 4, yA + 4) * B(xB + 4, yB + 4) + A(xA + 4, yA + 5) * B(xB + 5, yB + 4) + A(xA + 4, yA + 6) * B(xB + 6, yB + 4) + A(xA + 4, yA + 7) * B(xB + 7, yB + 4);
    C(xC + 4, yC + 5) += A(xA + 4, yA + 0) * B(xB + 0, yB + 5) + A(xA + 4, yA + 1) * B(xB + 1, yB + 5) + A(xA + 4, yA + 2) * B(xB + 2, yB + 5) + A(xA + 4, yA + 3) * B(xB + 3, yB + 5) + A(xA + 4, yA + 4) * B(xB + 4, yB + 5) + A(xA + 4, yA + 5) * B(xB + 5, yB + 5) + A(xA + 4, yA + 6) * B(xB + 6, yB + 5) + A(xA + 4, yA + 7) * B(xB + 7, yB + 5);
    C(xC + 4, yC + 6) += A(xA + 4, yA + 0) * B(xB + 0, yB + 6) + A(xA + 4, yA + 1) * B(xB + 1, yB + 6) + A(xA + 4, yA + 2) * B(xB + 2, yB + 6) + A(xA + 4, yA + 3) * B(xB + 3, yB + 6) + A(xA + 4, yA + 4) * B(xB + 4, yB + 6) + A(xA + 4, yA + 5) * B(xB + 5, yB + 6) + A(xA + 4, yA + 6) * B(xB + 6, yB + 6) + A(xA + 4, yA + 7) * B(xB + 7, yB + 6);
    C(xC + 4, yC + 7) += A(xA + 4, yA + 0) * B(xB + 0, yB + 7) + A(xA + 4, yA + 1) * B(xB + 1, yB + 7) + A(xA + 4, yA + 2) * B(xB + 2, yB + 7) + A(xA + 4, yA + 3) * B(xB + 3, yB + 7) + A(xA + 4, yA + 4) * B(xB + 4, yB + 7) + A(xA + 4, yA + 5) * B(xB + 5, yB + 7) + A(xA + 4, yA + 6) * B(xB + 6, yB + 7) + A(xA + 4, yA + 7) * B(xB + 7, yB + 7);
    C(xC + 5, yC + 0) += A(xA + 5, yA + 0) * B(xB + 0, yB + 0) + A(xA + 5, yA + 1) * B(xB + 1, yB + 0) + A(xA + 5, yA + 2) * B(xB + 2, yB + 0) + A(xA + 5, yA + 3) * B(xB + 3, yB + 0) + A(xA + 5, yA + 4) * B(xB + 4, yB + 0) + A(xA + 5, yA + 5) * B(xB + 5, yB + 0) + A(xA + 5, yA + 6) * B(xB + 6, yB + 0) + A(xA + 5, yA + 7) * B(xB + 7, yB + 0);
    C(xC + 5, yC + 1) += A(xA + 5, yA + 0) * B(xB + 0, yB + 1) + A(xA + 5, yA + 1) * B(xB + 1, yB + 1) + A(xA + 5, yA + 2) * B(xB + 2, yB + 1) + A(xA + 5, yA + 3) * B(xB + 3, yB + 1) + A(xA + 5, yA + 4) * B(xB + 4, yB + 1) + A(xA + 5, yA + 5) * B(xB + 5, yB + 1) + A(xA + 5, yA + 6) * B(xB + 6, yB + 1) + A(xA + 5, yA + 7) * B(xB + 7, yB + 1);
    C(xC + 5, yC + 2) += A(xA + 5, yA + 0) * B(xB + 0, yB + 2) + A(xA + 5, yA + 1) * B(xB + 1, yB + 2) + A(xA + 5, yA + 2) * B(xB + 2, yB + 2) + A(xA + 5, yA + 3) * B(xB + 3, yB + 2) + A(xA + 5, yA + 4) * B(xB + 4, yB + 2) + A(xA + 5, yA + 5) * B(xB + 5, yB + 2) + A(xA + 5, yA + 6) * B(xB + 6, yB + 2) + A(xA + 5, yA + 7) * B(xB + 7, yB + 2);
    C(xC + 5, yC + 3) += A(xA + 5, yA + 0) * B(xB + 0, yB + 3) + A(xA + 5, yA + 1) * B(xB + 1, yB + 3) + A(xA + 5, yA + 2) * B(xB + 2, yB + 3) + A(xA + 5, yA + 3) * B(xB + 3, yB + 3) + A(xA + 5, yA + 4) * B(xB + 4, yB + 3) + A(xA + 5, yA + 5) * B(xB + 5, yB + 3) + A(xA + 5, yA + 6) * B(xB + 6, yB + 3) + A(xA + 5, yA + 7) * B(xB + 7, yB + 3);
    C(xC + 5, yC + 4) += A(xA + 5, yA + 0) * B(xB + 0, yB + 4) + A(xA + 5, yA + 1) * B(xB + 1, yB + 4) + A(xA + 5, yA + 2) * B(xB + 2, yB + 4) + A(xA + 5, yA + 3) * B(xB + 3, yB + 4) + A(xA + 5, yA + 4) * B(xB + 4, yB + 4) + A(xA + 5, yA + 5) * B(xB + 5, yB + 4) + A(xA + 5, yA + 6) * B(xB + 6, yB + 4) + A(xA + 5, yA + 7) * B(xB + 7, yB + 4);
    C(xC + 5, yC + 5) += A(xA + 5, yA + 0) * B(xB + 0, yB + 5) + A(xA + 5, yA + 1) * B(xB + 1, yB + 5) + A(xA + 5, yA + 2) * B(xB + 2, yB + 5) + A(xA + 5, yA + 3) * B(xB + 3, yB + 5) + A(xA + 5, yA + 4) * B(xB + 4, yB + 5) + A(xA + 5, yA + 5) * B(xB + 5, yB + 5) + A(xA + 5, yA + 6) * B(xB + 6, yB + 5) + A(xA + 5, yA + 7) * B(xB + 7, yB + 5);
    C(xC + 5, yC + 6) += A(xA + 5, yA + 0) * B(xB + 0, yB + 6) + A(xA + 5, yA + 1) * B(xB + 1, yB + 6) + A(xA + 5, yA + 2) * B(xB + 2, yB + 6) + A(xA + 5, yA + 3) * B(xB + 3, yB + 6) + A(xA + 5, yA + 4) * B(xB + 4, yB + 6) + A(xA + 5, yA + 5) * B(xB + 5, yB + 6) + A(xA + 5, yA + 6) * B(xB + 6, yB + 6) + A(xA + 5, yA + 7) * B(xB + 7, yB + 6);
    C(xC + 5, yC + 7) += A(xA + 5, yA + 0) * B(xB + 0, yB + 7) + A(xA + 5, yA + 1) * B(xB + 1, yB + 7) + A(xA + 5, yA + 2) * B(xB + 2, yB + 7) + A(xA + 5, yA + 3) * B(xB + 3, yB + 7) + A(xA + 5, yA + 4) * B(xB + 4, yB + 7) + A(xA + 5, yA + 5) * B(xB + 5, yB + 7) + A(xA + 5, yA + 6) * B(xB + 6, yB + 7) + A(xA + 5, yA + 7) * B(xB + 7, yB + 7);
    C(xC + 6, yC + 0) += A(xA + 6, yA + 0) * B(xB + 0, yB + 0) + A(xA + 6, yA + 1) * B(xB + 1, yB + 0) + A(xA + 6, yA + 2) * B(xB + 2, yB + 0) + A(xA + 6, yA + 3) * B(xB + 3, yB + 0) + A(xA + 6, yA + 4) * B(xB + 4, yB + 0) + A(xA + 6, yA + 5) * B(xB + 5, yB + 0) + A(xA + 6, yA + 6) * B(xB + 6, yB + 0) + A(xA + 6, yA + 7) * B(xB + 7, yB + 0);
    C(xC + 6, yC + 1) += A(xA + 6, yA + 0) * B(xB + 0, yB + 1) + A(xA + 6, yA + 1) * B(xB + 1, yB + 1) + A(xA + 6, yA + 2) * B(xB + 2, yB + 1) + A(xA + 6, yA + 3) * B(xB + 3, yB + 1) + A(xA + 6, yA + 4) * B(xB + 4, yB + 1) + A(xA + 6, yA + 5) * B(xB + 5, yB + 1) + A(xA + 6, yA + 6) * B(xB + 6, yB + 1) + A(xA + 6, yA + 7) * B(xB + 7, yB + 1);
    C(xC + 6, yC + 2) += A(xA + 6, yA + 0) * B(xB + 0, yB + 2) + A(xA + 6, yA + 1) * B(xB + 1, yB + 2) + A(xA + 6, yA + 2) * B(xB + 2, yB + 2) + A(xA + 6, yA + 3) * B(xB + 3, yB + 2) + A(xA + 6, yA + 4) * B(xB + 4, yB + 2) + A(xA + 6, yA + 5) * B(xB + 5, yB + 2) + A(xA + 6, yA + 6) * B(xB + 6, yB + 2) + A(xA + 6, yA + 7) * B(xB + 7, yB + 2);
    C(xC + 6, yC + 3) += A(xA + 6, yA + 0) * B(xB + 0, yB + 3) + A(xA + 6, yA + 1) * B(xB + 1, yB + 3) + A(xA + 6, yA + 2) * B(xB + 2, yB + 3) + A(xA + 6, yA + 3) * B(xB + 3, yB + 3) + A(xA + 6, yA + 4) * B(xB + 4, yB + 3) + A(xA + 6, yA + 5) * B(xB + 5, yB + 3) + A(xA + 6, yA + 6) * B(xB + 6, yB + 3) + A(xA + 6, yA + 7) * B(xB + 7, yB + 3);
    C(xC + 6, yC + 4) += A(xA + 6, yA + 0) * B(xB + 0, yB + 4) + A(xA + 6, yA + 1) * B(xB + 1, yB + 4) + A(xA + 6, yA + 2) * B(xB + 2, yB + 4) + A(xA + 6, yA + 3) * B(xB + 3, yB + 4) + A(xA + 6, yA + 4) * B(xB + 4, yB + 4) + A(xA + 6, yA + 5) * B(xB + 5, yB + 4) + A(xA + 6, yA + 6) * B(xB + 6, yB + 4) + A(xA + 6, yA + 7) * B(xB + 7, yB + 4);
    C(xC + 6, yC + 5) += A(xA + 6, yA + 0) * B(xB + 0, yB + 5) + A(xA + 6, yA + 1) * B(xB + 1, yB + 5) + A(xA + 6, yA + 2) * B(xB + 2, yB + 5) + A(xA + 6, yA + 3) * B(xB + 3, yB + 5) + A(xA + 6, yA + 4) * B(xB + 4, yB + 5) + A(xA + 6, yA + 5) * B(xB + 5, yB + 5) + A(xA + 6, yA + 6) * B(xB + 6, yB + 5) + A(xA + 6, yA + 7) * B(xB + 7, yB + 5);
    C(xC + 6, yC + 6) += A(xA + 6, yA + 0) * B(xB + 0, yB + 6) + A(xA + 6, yA + 1) * B(xB + 1, yB + 6) + A(xA + 6, yA + 2) * B(xB + 2, yB + 6) + A(xA + 6, yA + 3) * B(xB + 3, yB + 6) + A(xA + 6, yA + 4) * B(xB + 4, yB + 6) + A(xA + 6, yA + 5) * B(xB + 5, yB + 6) + A(xA + 6, yA + 6) * B(xB + 6, yB + 6) + A(xA + 6, yA + 7) * B(xB + 7, yB + 6);
    C(xC + 6, yC + 7) += A(xA + 6, yA + 0) * B(xB + 0, yB + 7) + A(xA + 6, yA + 1) * B(xB + 1, yB + 7) + A(xA + 6, yA + 2) * B(xB + 2, yB + 7) + A(xA + 6, yA + 3) * B(xB + 3, yB + 7) + A(xA + 6, yA + 4) * B(xB + 4, yB + 7) + A(xA + 6, yA + 5) * B(xB + 5, yB + 7) + A(xA + 6, yA + 6) * B(xB + 6, yB + 7) + A(xA + 6, yA + 7) * B(xB + 7, yB + 7);
    C(xC + 7, yC + 0) += A(xA + 7, yA + 0) * B(xB + 0, yB + 0) + A(xA + 7, yA + 1) * B(xB + 1, yB + 0) + A(xA + 7, yA + 2) * B(xB + 2, yB + 0) + A(xA + 7, yA + 3) * B(xB + 3, yB + 0) + A(xA + 7, yA + 4) * B(xB + 4, yB + 0) + A(xA + 7, yA + 5) * B(xB + 5, yB + 0) + A(xA + 7, yA + 6) * B(xB + 6, yB + 0) + A(xA + 7, yA + 7) * B(xB + 7, yB + 0);
    C(xC + 7, yC + 1) += A(xA + 7, yA + 0) * B(xB + 0, yB + 1) + A(xA + 7, yA + 1) * B(xB + 1, yB + 1) + A(xA + 7, yA + 2) * B(xB + 2, yB + 1) + A(xA + 7, yA + 3) * B(xB + 3, yB + 1) + A(xA + 7, yA + 4) * B(xB + 4, yB + 1) + A(xA + 7, yA + 5) * B(xB + 5, yB + 1) + A(xA + 7, yA + 6) * B(xB + 6, yB + 1) + A(xA + 7, yA + 7) * B(xB + 7, yB + 1);
    C(xC + 7, yC + 2) += A(xA + 7, yA + 0) * B(xB + 0, yB + 2) + A(xA + 7, yA + 1) * B(xB + 1, yB + 2) + A(xA + 7, yA + 2) * B(xB + 2, yB + 2) + A(xA + 7, yA + 3) * B(xB + 3, yB + 2) + A(xA + 7, yA + 4) * B(xB + 4, yB + 2) + A(xA + 7, yA + 5) * B(xB + 5, yB + 2) + A(xA + 7, yA + 6) * B(xB + 6, yB + 2) + A(xA + 7, yA + 7) * B(xB + 7, yB + 2);
    C(xC + 7, yC + 3) += A(xA + 7, yA + 0) * B(xB + 0, yB + 3) + A(xA + 7, yA + 1) * B(xB + 1, yB + 3) + A(xA + 7, yA + 2) * B(xB + 2, yB + 3) + A(xA + 7, yA + 3) * B(xB + 3, yB + 3) + A(xA + 7, yA + 4) * B(xB + 4, yB + 3) + A(xA + 7, yA + 5) * B(xB + 5, yB + 3) + A(xA + 7, yA + 6) * B(xB + 6, yB + 3) + A(xA + 7, yA + 7) * B(xB + 7, yB + 3);
    C(xC + 7, yC + 4) += A(xA + 7, yA + 0) * B(xB + 0, yB + 4) + A(xA + 7, yA + 1) * B(xB + 1, yB + 4) + A(xA + 7, yA + 2) * B(xB + 2, yB + 4) + A(xA + 7, yA + 3) * B(xB + 3, yB + 4) + A(xA + 7, yA + 4) * B(xB + 4, yB + 4) + A(xA + 7, yA + 5) * B(xB + 5, yB + 4) + A(xA + 7, yA + 6) * B(xB + 6, yB + 4) + A(xA + 7, yA + 7) * B(xB + 7, yB + 4);
    C(xC + 7, yC + 5) += A(xA + 7, yA + 0) * B(xB + 0, yB + 5) + A(xA + 7, yA + 1) * B(xB + 1, yB + 5) + A(xA + 7, yA + 2) * B(xB + 2, yB + 5) + A(xA + 7, yA + 3) * B(xB + 3, yB + 5) + A(xA + 7, yA + 4) * B(xB + 4, yB + 5) + A(xA + 7, yA + 5) * B(xB + 5, yB + 5) + A(xA + 7, yA + 6) * B(xB + 6, yB + 5) + A(xA + 7, yA + 7) * B(xB + 7, yB + 5);
    C(xC + 7, yC + 6) += A(xA + 7, yA + 0) * B(xB + 0, yB + 6) + A(xA + 7, yA + 1) * B(xB + 1, yB + 6) + A(xA + 7, yA + 2) * B(xB + 2, yB + 6) + A(xA + 7, yA + 3) * B(xB + 3, yB + 6) + A(xA + 7, yA + 4) * B(xB + 4, yB + 6) + A(xA + 7, yA + 5) * B(xB + 5, yB + 6) + A(xA + 7, yA + 6) * B(xB + 6, yB + 6) + A(xA + 7, yA + 7) * B(xB + 7, yB + 6);
    C(xC + 7, yC + 7) += A(xA + 7, yA + 0) * B(xB + 0, yB + 7) + A(xA + 7, yA + 1) * B(xB + 1, yB + 7) + A(xA + 7, yA + 2) * B(xB + 2, yB + 7) + A(xA + 7, yA + 3) * B(xB + 3, yB + 7) + A(xA + 7, yA + 4) * B(xB + 4, yB + 7) + A(xA + 7, yA + 5) * B(xB + 5, yB + 7) + A(xA + 7, yA + 6) * B(xB + 6, yB + 7) + A(xA + 7, yA + 7) * B(xB + 7, yB + 7);
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
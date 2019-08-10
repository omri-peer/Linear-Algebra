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

    static Matrix add(const Matrix& A, const Matrix& B, int xA, int yA, int xB, int yB, int size)
    {
        Matrix output(size, size);
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                output.major_row[i][j] = A.major_row[xA + i][yA + j] + B.major_row[xB + i][yB + j];
            }
        }

        return std::move(output);
    }

    static Matrix sub(const Matrix& A, const Matrix& B, int xA, int yA, int xB, int yB, int size)
    {
        Matrix output(size, size);
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                output.major_row[i][j] = A.major_row[xA + i][yA + j] - B.major_row[xB + i][yB + j];
            }
        }

        return std::move(output);
    }

    static void add_in_place(Matrix& A, const Matrix& B, int xA, int yA, int xB, int yB, int size)
    {
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                A.major_row[xA + i][yA + j] += B.major_row[xB + i][yB + j];
            }
        }
    }

    static void sub_in_place(Matrix& A, const Matrix& B, int xA, int yA, int xB, int yB, int size)
    {
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                A.major_row[xA + i][yA + j] -= B.major_row[xB + i][yB + j];
            }
        }
    }

    static void copy(Matrix& A, const Matrix& B, int xA, int yA, int xB, int yB, int size)
    {
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                A.major_row[xA + i][yA + j] = B.major_row[xB + i][yB + j];
            }
        }
    }

    static void mul(const Matrix& A, const Matrix& B, Matrix& C, int xA, int yA, int xB, int yB, int xC, int yC, int size)
    {
        //        std::cout << "A1\n" << A << std::endl;
        if (size == 8) {
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
            return;
        }

        int half_size = size / 2;
        //        Matrix<T> A11(half_size, half_size);
        //        Matrix<T> A12(half_size, half_size);
        //        Matrix<T> A21(half_size, half_size);
        //        Matrix<T> A22(half_size, half_size);
        //        Matrix<T> B11(half_size, half_size);
        //        Matrix<T> B12(half_size, half_size);
        //        Matrix<T> B21(half_size, half_size);
        //        Matrix<T> B22(half_size, half_size);

        //        Matrix<T> M1(half_size, half_size);
        //        Matrix<T> M2(half_size, half_size);
        //        Matrix<T> M3(half_size, half_size);
        //        Matrix<T> M4(half_size, half_size);
        //        Matrix<T> M5(half_size, half_size);
        //        Matrix<T> M6(half_size, half_size);
        //        Matrix<T> M7(half_size, half_size);
        //        std::cout << "A2\n" << A << std::endl;
        //        for (int i = 0; i < size; ++i) {
        //            for (int j = 0; j < size; ++j) {
        //                if (i < half_size && j < half_size) {
        //                    //                    A11(i, j) = A(xA + i, yA + j);
        //                    B11(i, j) = B(xB + i, yB + j);
        //                }
        //                if (i < half_size && j >= half_size) {
        //                    A12(i, j - half_size) = A(xA + i, yA + j);
        //                    //                    B12(i, j - half_size) = B(xB + i, yB + j);
        //                }
        //                if (i >= half_size && j < half_size) {
        //                    A21(i - half_size, j) = A(xA + i, yA + j);
        //                    B21(i - half_size, j) = B(xB + i, yB + j);
        //                }
        //                if (i >= half_size && j >= half_size) {
        //                    A22(i - half_size, j - half_size) = A(xA + i, yA + j);
        //                    B22(i - half_size, j - half_size) = B(xB + i, yB + j);
        //                }
        //            }
        //        }
        // if (size >= 256) std::cout << size << std::endl;

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

        /*
                mul(A12 - A22, B21 + B22, C, 0, 0, 0, 0, xC, yC, half_size);
                B21 -= B11;
                mul(A, B21, C, xA + half_size, yA + half_size, 0, 0, xC + half_size, yC, half_size);
                mul(add(A22, A, 0, 0, xA, yA, half_size), B22 + B11, C, 0, 0, 0, 0, xC + half_size, yC + half_size, half_size);
                A22 += A21;
                Matrix M(half_size, half_size);
                mul(A22, B, M, 0, 0, xB, yB, 0, 0, half_size);
                A22 = M;
                sub_in_place(B22, B, 0, 0, xB, yB + half_size, half_size);
                mul(A, B22, M, xA, yA, 0, 0, 0, 0, half_size);
                B22 = M;
                add_in_place(A12, A, 0, 0, xA, yA, half_size);
                mul(A12, B, C, 0, 0, xB + half_size, yB + half_size, xC, yC + half_size, half_size);

                add_in_place(C, C, xC, yC, xC + half_size, yC + half_size, half_size);
                add_in_place(C, C, xC, yC, xC + half_size, yC, half_size);
                sub_in_place(C, C, xC, yC, xC, yC + half_size, half_size);
                sub_in_place(C, B22, xC, yC + half_size, 0, 0, half_size);
                add_in_place(C, A22, xC + half_size, yC, 0, 0, half_size);
                sub_in_place(A21, A, 0, 0, xA, yA, half_size);
                add_in_place(B11, B, 0, 0, xB, yB + half_size, half_size);
                mul(A21, B11, M, 0, 0, 0, 0, 0, 0, half_size);
                A21 = M;
                B22 -= A21;
                A22 += B22;
                sub_in_place(C, A22, xC + half_size, yC + half_size, 0, 0, half_size);
                */
        //        std::cout << "A3\n" << A << std::endl;
        /*
        mul(A12 - A22, B21 + B22, M6, 0, 0, 0, 0, 0, 0, half_size);
        B21 -= B11;
        mul(A, B21, M4, xA + half_size, yA + half_size, 0, 0, 0, 0, half_size);
        mul(A22 + A11, B22 + B11, M1, 0, 0, 0, 0, 0, 0, half_size);
        A22 += A21;
        mul(A22, B, M2, 0, 0, xB, yB, 0, 0, half_size);
        B22 -= B12;
        mul(A, B22, M3, xA, yA, 0, 0, 0, 0, half_size);
        A12 += A11;
        mul(A12, B22, M5, 0, 0, 0, 0, 0, 0, half_size);
        //        std::cout << "A4\n" << A << std::endl;
        M6 += M1;
        M6 += M4;
        M6 -= M5;
        M5 -= M3;
        M4 += M2;
        A21 -= A11;
        B11 += B12;
        mul(A21, B11, M7, 0, 0, 0, 0, 0, 0, half_size);
        M3 -= M7;
        M2 += M3;
        M1 -= M2;*/
        //        std::cout << "A5\n" << A << std::endl;
        /*
                M6 = (A12 - A22) * (B12 + B22);
                B21 -= B11;
                M4 = A22 * B21;
                M1 = (A22 + A11) * (B22 + B11);
                A22 += A21;
                M2 = A22 * B11;
                B22 -= B12;
                M3 = A11 * B22;
                A12 += A11;
                M5 = A12 * B22;

                M6 += M1;
                M6 += M4;
                M6 -= M5;
                M5 -= M3;
                M4 += M2;
                A21 -= A11;
                B11 += B12;
                M7 = A21 * B11;
                M3 -= M7;
                M2 += M3;
                M1 -= M2;
        */
        //        std::cout << "M6\n" << M6 << std::endl;
        //        std::cout << "M5\n" << M5 << std::endl;
        //        std::cout << "M4\n" << M4 << std::endl;
        //        std::cout << "M1\n" << M1 << std::endl;
        /*
                for (int i = 0; i < size; i++) {
                    for (int j = 0; j < size; j++) {
                        if (i < half_size && j < half_size) {
                            C(xC + i, yC + j) = M6(i, j);
                        }
                        if (i < half_size && j >= half_size) {
                            C(xC + i, yC + j) = M5(i, j - half_size);
                        }
                        if (i >= half_size && j < half_size) {
                            C(xC + i, yC + j) = M4(i - half_size, j);
                        }
                        if (i >= half_size && j >= half_size) {
                            C(xC + i, yC + j) = M1(i - half_size, j - half_size);
                        }
                    }
                }
        */
        /*
        Matrix A1(size, size);
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                A1(i, j) = A(xA + i, yA + j);
            }
        }

        Matrix B1(size, size);
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                B1(i, j) = B(xB + i, yB + j);
            }
        }

        Matrix C1(size, size);
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                C1(i, j) = C(xC + i, yC + j);
            }
        }

        if (A1 * B1 != C1) {
            std::cout << "!!!!!!!!!!!!!!\n" << A1 << std::endl << B1 << std::endl << C1 << std::endl;
        }
         */
        //        std::cout << size << std::endl;
    }

    // gets one square matrix
    Matrix<T> gaussian_elimination(T& det) const
    {
        int size = (*this).get_cols();
        Matrix<T> copy(*this);
        Matrix<T> inverse(size, size);
        // set inverse to unit matrix
        for (int i = 0; i < size; ++i) {
            inverse(i, i) = 1;
        }

        for (int row = 0; row < size; ++row) {
            // if copy(i,i) is 0, swap it with a lower row.
            for (int lower_row = row + 1; lower_row < size && copy(row, row) == 0; ++lower_row) {
                if (copy(lower_row, row) != 0) {
                    copy.swap_rows(row, lower_row);
                    inverse.swap_rows(row, lower_row);
                    det *= (-1);
                }
            }

            if (copy(row, row) != 0) {
                T scalar = copy(row, row);
                copy.multiply_row_by_scalar(row, (1 / scalar));
                inverse.multiply_row_by_scalar(row, (1 / scalar));
                det *= scalar;
                // set all other rows to zero in the row-th index
                for (int other_row = 0; other_row < size; ++other_row) {
                    if (other_row != row && copy(other_row, row) != 0) {
                        scalar = copy(other_row, row);
                        copy.add_multiplied_row(other_row, row, scalar * (-1));
                        inverse.add_multiplied_row(other_row, row, scalar * (-1));
                    }
                }
            }
            else {
                det = 0;
            }
        }
        return inverse;
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
                major_row[i][j] += m.major_row[i][j];
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
        return std::move(sum);
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
            return paddedProd;
        }

        if (size == 2) {
            paddedProd(0, 0) = A(0, 0) * B(0, 0) + A(0, 1) * B(1, 0);
            paddedProd(0, 1) = A(0, 0) * B(0, 1) + A(0, 1) * B(1, 1);
            paddedProd(1, 0) = A(1, 0) * B(0, 0) + A(1, 1) * B(1, 0);
            paddedProd(1, 1) = A(1, 0) * B(0, 1) + A(1, 1) * B(1, 1);
            return paddedProd;
        }

        if (size == 4) {
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
        return std::move(neg);
    }

    // in-place matrices substraction
    void operator-=(const Matrix& m)
    {
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; ++j) {
                major_row[i][j] -= m.major_row[i][j];
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
        return std::move(diff);
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
                major_row[i][j] *= s;
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
        return std::move(prod);
    }

    friend Matrix operator*(const T& s, const Matrix& m)
    {
        return m * s;
    }

    // gets two rows and swap them, the rows should have a valid range
    void swap_rows(int i, int j)
    {
        vector_t temp(std::move(major_row[i]));
        major_row[i] = std::move(major_row[j]);
        major_row[j] = std::move(temp);
    }

    void multiply_row_by_scalar(int row, T scalar)
    {
        for (int i = 0; i < (*this).get_cols(); ++i) {
            (*this)(row, i) *= scalar;
        }
    }

    void add_multiplied_row(int row1, int row2, T scalar)
    {
        for (int i = 0; i < (*this).get_cols(); ++i) {
            (*this)(row1, i) += (*this)(row2, i) * scalar;
        }
    }

    // gets one square matrix
    T det() const
    {
        T det = 1;
        (*this).gaussian_elimination(det);
        return det;
    }

    Matrix find_inverse() const
    {
        T det;
        return (*this).gaussian_elimination(det);
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
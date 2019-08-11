#pragma once

#include <iostream>
#include <vector>
template <class T>
class Vector {
    using vector_t = std::vector<T>;

private:
    vector_t entries; // entries store the vector's entries, 0 based

public:
    // default constructor. Constructs a vector of size 1 with 0
    explicit Vector() : entries(vector_t(1, 0))
    {
    }

    // parametrized constructor. Constructs a vector of given size, with 0's
    // assumes size is positive
    explicit Vector(unsigned int size) : entries(size, 0)
    {
    }

    // constructor by const reference to std::vector of entries
    explicit Vector(const vector_t& ents) : entries(ents)
    {
    }

    // constructor by rvalue reference to std::vector of entries
    explicit Vector(vector_t&& ents) : entries(std::move(ents))
    {
    }

    // copy constructor by const reference
    Vector(const Vector& v) : entries(v.entries)
    {
    }

    // move constructor by rvalue reference to std::vector of entries
    Vector(Vector&& v) noexcept : entries(std::move(v.entries))
    {
    }

    // copy assignment by const reference
    Vector& operator=(const Vector& v)
    {
        entries = v.entries;
    }

    // move assignment by rvalue reference
    Vector& operator=(Vector&& v) noexcept
    {
        entries = std::move(v.entries);
    }

    // returns whether the vectors have the same entries or not
    bool operator==(const Vector& v) const
    {
        return entries == v.entries;
    }
    // returns whether the vectors have different entries or not
    bool operator!=(const Vector& v) const
    {
        return !(*this == v);
    }

    // returns the size of a Vector (number of entries)
    unsigned int dimension() const
    {
        return entries.size();
    }

    // v(i) is the i-th entry of the vector (for get and set)
    // both the getter and the setter assume the given index i satisfies 0 <= i < (*this).dimension()
    T operator()(unsigned int i) const
    {
        return entries[i];
    }
    T& operator()(unsigned int i)
    {
        return entries[i];
    }

    // in-place multiplication by a scalar
    void operator*=(T scalar)
    {
        for (unsigned int i = 0; i < dimension(); i++)
            (*this)(i) *= scalar;
    }

    // multiplication by a scalar
    friend Vector operator*(const Vector& v, T scalar)
    {
        Vector new_vector(v.entries);
        for (unsigned int i = 0; i < v.dimension(); i++)
            new_vector(i) *= scalar;

        return std::move(new_vector);
    }

    // reversed multiplication by a scalar (exactly the same because vector by scalar multiplication is a commutative binary operation)
    friend Vector operator*(T scalar, const Vector& v)
    {
        return v * scalar;
    }

    // inner multiplication
    // assumes other's dimension is the same as this vector's
    T operator*(const Vector& other) const
    {
        T res = 0;
        for (unsigned int i = 0; i < dimension(); i++)
            res += (*this)(i)*other(i);

        return res;
    }

    // in-place vector addition
    // assumes other's dimension is the same as this vector's
    void operator+=(const Vector& other)
    {
        for (unsigned int i = 0; i < dimension(); i++)
            (*this)(i) += other(i);
    }

    // vector addition
    // assumes other's dimension is the same as this vector's
    Vector operator+(const Vector& other) const
    {
        Vector new_vector(entries);
        for (unsigned int i = 0; i < dimension(); i++)
            new_vector(i) += other(i);

        return std::move(new_vector);
    }

    // vector negation
    Vector operator-() const
    {
        Vector new_vector(entries);
        for (unsigned int i = 0; i < dimension(); i++)
            new_vector(i) = -new_vector(i);

        return std::move(new_vector);
    }

    // in-place vector subtraction
    // assumes other's dimension is the same as this vector's
    void operator-=(const Vector& other)
    {
        for (unsigned int i = 0; i < dimension(); i++)
            (*this)(i) -= other(i);
    }

    // vector subtraction
    // assumes other's dimension is the same as this vector's
    Vector operator-(const Vector& other) const
    {
        Vector new_vector(entries);
        for (unsigned int i = 0; i < dimension(); i++)
            new_vector(i) -= other(i);

        return std::move(new_vector);
    }
};

// sending to output stream using <<
template <class T>
std::ostream& operator<<(std::ostream& strm, const Vector<T>& v)
{
    strm << "(";
    for (unsigned int i = 0; i < v.dimension() - 1; i++) {
        strm << v(i) << ", ";
    }
    strm << v(v.dimension() - 1) << ")\n";

    return strm;
}

// apply the gram - schmidt algorithm to a given basis - in order to return a new, orthogonal basis, every prefix of which spans the same subspace as the corresponding prefix of the given basis
// implemented in the naive way
// assumes all vectors in the given basis are of the same dimension, and that they are all linearly independent
template <class T>
std::vector<Vector<T>> gs(const std::vector<Vector<T>>& basis)
{
    std::vector<Vector<T>> new_basis;
    for (unsigned int i = 0; i < basis.size(); i++) {
        new_basis.push_back(basis[i]);
        for (unsigned int j = 0; j < i; j++) {
            new_basis[i] -= Vector<T>(new_basis[j] * (basis[i] * new_basis[j] / (new_basis[j] * new_basis[j])));
        }
    }

    return std::move(new_basis);
}
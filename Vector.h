#pragma once

#include <iostream>
#include <vector>
template <class T>
class Vector {
    using vector_t = std::vector<T>; // entries store the vector's entries, 0 based

private:
    vector_t entries;

public:
    // default constructor. Constructs a vector of size 1 with 0
    explicit Vector() : entries(vector_t(1, 0))
    {
    }

    // parametrized constructor. Constructs a vector of given size, with 0's
    explicit Vector(int size) : entries(size, 0)
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

    bool operator!=(const Vector& v) const
    {
        return !(*this == v);
    }

    // returns the size of a Vector (number of entries)
    int dimension() const
    {
        return entries.size();
    }

    // v(i) is the i-th entry of the vector (for get and set)
    T operator()(int i) const
    {
        return entries[i];
    }
    T& operator()(int i)
    {
        return entries[i];
    }

    // in-place multiplication by a scalar
    void operator*=(T scalar)
    {
        for (int i = 0; i < dimension(); i++)
            (*this)(i) *= scalar;
    }

    // multiplication by a scalar
    friend Vector operator*(const Vector& v, T scalar)
    {
        Vector new_vector(v.entries);
        for (int i = 0; i < v.dimension(); i++)
            new_vector(i) *= scalar;

        return new_vector;
    }

    // reversed multiplication by a scalar
    friend Vector operator*(T scalar, const Vector& v)
    {
        return v * scalar;
    }

    // inner multiplication
    T operator*(const Vector& other) const
    {
        T res = 0;
        for (int i = 0; i < dimension(); i++)
            res += (*this)(i)*other(i);

        return res;
    }

    // in-place vectors addition
    void operator+=(const Vector& other)
    {
        for (int i = 0; i < dimension(); i++)
            (*this)(i) += other(i);
    }

    // vectors addition
    Vector operator+(const Vector& other) const
    {
        Vector new_vector(entries);
        for (int i = 0; i < dimension(); i++)
            new_vector(i) += other(i);

        return new_vector;
    }

    // vector negation
    Vector operator-() const
    {
        Vector new_vector(entries);
        for (int i = 0; i < dimension(); i++)
            new_vector(i) = -new_vector(i);

        return new_vector;
    }

    // in-place vectors subtraction
    void operator-=(const Vector& other)
    {
        for (int i = 0; i < dimension(); i++)
            (*this)(i) -= other(i);
    }

    // vectors subtraction
    Vector operator-(const Vector& other) const
    {
        Vector new_vector(entries);
        for (int i = 0; i < dimension(); i++)
            new_vector(i) -= other(i);

        return new_vector;
    }
};

// sending to output stream using <<
template <class T>
std::ostream& operator<<(std::ostream& strm, const Vector<T>& v)
{
    strm << "(";
    for (int i = 0; i < v.dimension() - 1; i++) {
        strm << v(i) << ", ";
    }
    strm << v(v.dimension() - 1) << ")\n";

    return strm;
}

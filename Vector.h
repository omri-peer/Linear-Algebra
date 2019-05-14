#pragma once

#include <vector>
#include <iostream>
template <class T>
class Vector
{
    using vector_t = std::vector<T>;

private:
    vector_t entries;

public:
    Vector<T>()
    {
        entries = vector_t(1, 0);
    }

    Vector<T>(int size) : entries(size, 0)
    {
    }

    Vector<T>(const vector_t& v)
    {
        entries = std::move(v);
    }

    Vector<T>(const Vector& v)
    {
        entries = std::move(v.entries);
    }

    int dimension() const
    {
        return entries.size();
    }

		T operator()(int  i) const
		{
			return entries[i];
		}
		T& operator()(int  i)
		{
			return entries[i];
		}

		template <class U>
		friend Vector<U> operator*(const Vector<U>& v, U scalar);

		template <class U>
		friend Vector<U> operator*(U scalar, const Vector<U>& v);

		T operator*(const Vector<T>& other) const
		{
			T res = 0;
			for (int i = 0; i < dimension(); i++)
				res += entries[i] * other(i);

			return res;
		}

		Vector<T> operator+(const Vector& other) const
		{
            vector_t new_entries = entries;
            for (int i = 0; i < dimension(); i++)
                new_entries[i] += other(i);

            return new_entries;
        }

        Vector<T> operator-() const
		{
            vector_t new_entries = entries;
            for (int i = 0; i < dimension(); i++)
                new_entries[i] = -new_entries[i];

            return new_entries;
        }

        Vector<T> operator-(const Vector& other) const
		{
            vector_t new_entries = entries;
            for (int i = 0; i < dimension(); i++)
                new_entries[i] -= other(i);

            return new_entries;
        }
};

template <class T>
Vector<T> operator*(const Vector<T>& v, T scalar)
{
    std::vector<T> entries = v.entries;
    for (int i = 0; i < v.dimension(); i++)
        entries[i] *= scalar;

    return entries;
}

template <class T>
Vector<T> operator*(T scalar, const Vector<T>& v)
{
    std::vector<T> entries = v.entries;
    for (int i = 0; i < v.dimension(); i++)
        entries[i] *= scalar;

    return entries;
}

template <class T>
std::ostream& operator<<(std::ostream& strm, const Vector<T>& v)
{
	strm << "(";
	for (int i = 0; i < v.dimension()-1; i++)
	{
		strm << v(i) << ", ";
	}
	strm << v(v.dimension()-1) << ")\n";

	return strm;
}

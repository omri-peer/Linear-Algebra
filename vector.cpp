#include <iostream>
#include "Matrix.h"

int main()
{
	/*std::vector<double> ents;
	ents.push_back(2);
	ents.push_back(0);
	ents.push_back(-1);

	Vector<double> v1 = Vector<double>(ents);
	Vector<double> v2 = Vector<double>(std::vector<double>(3, 2));

	std::cout << v1;
	std::cout << v2;
	std::cout << v1*4.0;
	std::cout << (-0.5)*v2;
	std::cout << v1 + 2.0*v2;
	std::cout << -v2;
	std::cout << v1*v2 << "\n";
	std::cout << (v1 - v2);

	v1 = v2;
	v2 = Vector<double>(ents);

	std::cout << v1;
	std::cout << v2;*/
	std::vector<double> ents;
	ents.push_back(2);
	ents.push_back(0);
	ents.push_back(-1);

	Vector<double> v1 = Vector<double>(ents);

	Matrix<double> m1(2,2);
	m1(1, 1) = 4;
	m1(0, 0) = 1;
	m1(0, 1) = 2;
	Matrix<double> m2=m1+m1;

	Matrix<double> m3=m1*m2;

	std::cout << m1;
	std::cout << m2;
	std::cout << m3;
}
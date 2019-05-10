#include <iostream>
#include <vector>
#include <math.h>
#include <cmath>
#include "Matrix.h"
#include "Vector.h"

std::vector<Vector<double> > gs(const std::vector<Vector<double> > &basis)
{
	std::vector<Vector<double> > new_basis;
	for (int i = 0; i < basis[0].dimension(); i++)
	{
		new_basis.push_back(basis[i]);
		for (int j = 0; j < i; j++)
		{
			new_basis[i] = new_basis[i] - Vector<double>(new_basis[j]*(basis[i]*new_basis[j]/(new_basis[j]*new_basis[j])));
		}
	}

	return new_basis;
}

std::vector<Vector<double> > update(int size, Matrix<double> & m, const std::vector<Vector<double> > &basis)
{
	std::vector<Vector<double> > ortho = gs(basis);

	for (int i=0; i<size; i++)
	{
		for (int j=0; j<size; j++)
		{
			m(i,j)=(basis[i]*ortho[j])/(ortho[j]*ortho[j]);
		}
	}
	
	return ortho;
}

std::vector<Vector<double> > lll(double delta, const std::vector<Vector<double> > &basis)
{
	Vector<double> temp;

	std::vector<Vector<double> > our_basis = basis;
	int size = basis.size();

	Matrix<double> m(size, size);
	std::vector<Vector<double> > ortho = update(size, m, our_basis);

	int k = 1;
	while(k <= size - 1)
	{
		for (int j = k - 1; j >= 0; --j)
		{
			if (std::abs(m(k, j)) > 0.5)
			{
				our_basis[k] = our_basis[k] - floor(0.5 + m(k, j)) * our_basis[j];
				ortho = update(size, m, our_basis);
			}
		}
		if (ortho[k] * ortho[k] >= (delta - m(k, k - 1) * m(k, k - 1)) * (ortho[k-1] * ortho[k-1]))
		{
			++k;
		}
		else
		{
			temp = our_basis[k];
			our_basis[k] = our_basis[k - 1];
			our_basis[k - 1] = temp;

			ortho = update(size, m, our_basis);

			if (k - 1 > 1) --k;
			else k = 1;
		}
	}

	return our_basis;
}


int main()
{
	std::vector<double> ent1(3, 0);
	ent1[0] = 0;
	ent1[1] = 10;
	ent1[2] = 2;
	std::vector<double> ent2(3, 200);
	ent2[1] = 314;
	ent2[2] = 503;
	std::vector<double> ent3(3, 20);
	ent3[1] = 30;
	ent3[2] = 50;

	Vector<double> v1(ent1);
	Vector<double> v2(ent2);
	Vector<double> v3(ent3);

	std::vector<Vector<double> > basis;
	basis.push_back(v1);
	basis.push_back(v2);
	basis.push_back(v3);

	std::cout << basis[0] << "\n" << basis[1] << "\n" << basis[2] << "\n";
	std::cout << "\n";

	std::vector<Vector<double> > lll_basis = lll(0.75, basis);
	std::cout << lll_basis[0] << "\n" << lll_basis[1] << "\n" << lll_basis[2] << "\n";

	return 0;
}
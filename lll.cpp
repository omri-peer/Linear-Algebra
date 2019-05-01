#include <iostream>
#include "Matrix.h"

double abs(double r)
{
	if (r>0)
	{
		return r;
	}
	return -r;
}

double round(double r)
{
	int inted= (int) r;
	if (r>inted)
	{
		if (r>=inted+0.5)
		{
			return inted+1;
		}
		return inted;
	}
	else
	{
		if (r>=inted-0.5)
		{
			return inted;
		}
		return inted-1;
	}
}

void gs(const std::vector<Vector<double> > &basis, std::vector<Vector<double> > &new_basis)
{
	for (int i = 0; i < basis[0].dimension(); i++)
	{
		new_basis.push_back(basis[i]);
		for (int j = 0; j < i; j++)
		{
			new_basis[i] = new_basis[i] - Vector<double>(basis[j]*(basis[i]*basis[j]/(basis[j]*basis[j])));
			std::cout << i << " " << j << " " << new_basis[i];
		}
	}

	return;
}

void update(int size, Matrix<double> & m, const std::vector<Vector<double> > &basis, std::vector<Vector<double> > &ortho)
{
	gs(basis, ortho);
	for (int i=0; i<size; i++)
	{
		for (int j=0; j<size; j++)
		{
			m(i,j)=(basis[i]*ortho[j])/(ortho[j]*ortho[j]);
		}
	}
	return;
}

std::vector<Vector<double> > lll(double delta, const std::vector<Vector<double> > &basis)
{
	Vector<double> temp;

	std::vector<Vector<double> > our_basis = basis;
	int size = basis.size();
	Matrix<double> m(size, size);
	std::vector<Vector<double> > ortho;
	update(size, m, our_basis, ortho);
	int k = 1;
	while(k <= size - 1)
	{
		for (int j = k - 1; j >= 0; --j)
		{
			if (abs(m(k, j)) > 0.5)
			{
				our_basis[k] = our_basis[k] - round(m(k, j)) * our_basis[j];
				update(size, m, our_basis, ortho);
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
			update(size, m, our_basis, ortho);
			if (k - 1 > 1) --k;
			else k = 1;
		}
	}
	return our_basis;
}


int main()
{
	std::vector<double> ent1(3, 0);
	ent1[0] = 1;
	std::vector<double> ent2(3, 0);
	ent2[1] = 1;
	std::vector<double> ent3(3, 1);

	Vector<double> v1(ent1);
	Vector<double> v2(ent2);
	Vector<double> v3(ent3);

	std::vector<Vector<double> > basis;
	basis.push_back(v1);
	basis.push_back(v2);
	basis.push_back(v3);

	std::vector<Vector<double> > new_basis;
	gs(basis, new_basis);
	std::cout << new_basis[0] << new_basis[1] << new_basis[2];
	return 0;
}
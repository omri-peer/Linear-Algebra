#pragma once

#include "Matrix.h"
#include "Vector.h"
#include <cmath>
#include <gmp.h>
#include <gmpxx.h>
#include <iostream>
#include <vector>

std::vector<Vector<mpq_class>> lll(double delta, const std::vector<Vector<mpq_class>>& basis);

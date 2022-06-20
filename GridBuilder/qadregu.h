#pragma once
#include <cmath>
static struct QuadrEqu
{
	static const double eps;
	static unsigned char solve(double a, double b, double c, double& x1, double& x2)
	{
		double descr = b * b - 4 * a * c;

		if (abs(descr / (b * b)) < eps)
		{
			x1 = -b / (2 * a);
		}

		if (descr < 0)
			return 0;

		descr = sqrt(descr);
		x1 = (-b - descr) / (2 * a);
		x2 = (-b + descr) / (2 * a);
		return 2;
	}
};
const double QuadrEqu::eps = 1e-14;

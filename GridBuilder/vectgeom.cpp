#include "vectgeom.h"
#include <cmath>

Vect multVectVect(Vect vec1, Vect vec2)
{
	Vect res;
	res.x = vec1.y * vec2.z - vec1.z * vec2.y;
	res.y = vec1.z * vec2.x - vec1.x * vec2.z;
	res.z = vec1.x * vec2.y - vec1.y * vec2.x;

	return res;
}

Vect normalizeVect(Vect vec)
{
	return vec / sqrt(multScalVect(vec, vec));
}

double multScalVect(Vect vec1, Vect vec2)
{
	Vect mult = vec1 * vec2;
	return mult.x + mult.y + mult.z;
}

Vect findNormVect(Coord p1, Coord p2, Coord p3)
{
	Vect vec1 = p2 - p1;
	Vect vec2 = p3 - p1;

	Vect vectMult = multVectVect(vec1, vec2);
	vectMult = normalizeVect(vectMult);

	return vectMult;
}
#pragma once
#include "coord.h"

typedef Coord Vect;

Vect multVectVect(Vect vec1, Vect vec2);
Vect normalizeVect(Vect vec);
double multScalVect(Vect vec1, Vect vec2);
Coord findNormVect(Coord p1, Coord p2, Coord p3);

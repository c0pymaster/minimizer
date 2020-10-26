#ifndef MELZAK_H
#define MELZAK_H

#include "interval.hpp"
#include "imath.hpp"
#include "geometry.h"
#include <math.h>
#include <assert.h>

using namespace std;
using namespace cxsc;

// This is the implementation of the Melzak algorithm with additional assumption that the number of given points is 4.

// Note that the reconstruction phase is considered successful if the result of some checks could not be determined
// (that usually happens if some points are very close to each other).
// It means that the algorithm can return the length smaller that the real length of Steiner tree on given points.
// This is fine since we are only interested in the lower bound on this length.

// Steiner tree on 4 points contains no more than 6 points (4 given points + 2 Steiner points)
const int maxn = 6;

struct Melzak {
	point p[maxn];

	// the reconstruction step of the Melzak algorithm;
	// returns true if the reconstruction is successfull
	// (also returns true if the result of some of the checks could not be determined)
	bool reconstruct(int s, int a, int b, int v) {
		// construct the circle through points p[a], p[b], p[s]
		point o;
		interval r;
		getCircle(p[a], p[b], p[s], o, r);

		if (getSign((p[s] - p[v]).len()) == 0) {
			// p[s] and p[v] are too close to each other to construct a line through them;
			// assume that the reconstruction is successful
			return true;
		}

		// find the intersection of the circle and the line which goes through p[s] and p[v]
		point res[2];
		intersect(line(p[s], p[v]), o, r, res);

		line l(p[a], p[b]);

		// if one of the intersection points is inside both the segment between p[s] and p[v], and the smaller arc between p[a] and p[b],
		// then this point is the desired Steiner point; otherwise, the reconstruction phase failed
		bool success = false;
		for (int i = 0; i < 2; ++i) {
			// if (p[s] - res[i]) * (p[v] - res[i]) > 0, then res[i] is not inside the segment between p[s] and p[v];
			// if res[i] and p[s] are on the same side of the line l, then res[i] is not inside the smaller arc between p[a] and p[b]
			if (getSign((p[s] - res[i]) * (p[v] - res[i])) <= 0 && l.getSide(p[s]) * l.getSide(res[i]) != 1) {
				success = true;

				// the total length of the tree should not change during the reconstruction phase
				assert(getSign((p[s] - res[i]).len() - (p[a] - res[i]).len() - (p[b] - res[i]).len()) == 0);

				p[s] = res[i];
				break;
			}
		}

		return success;
	}

	// return the third vertex of the equilateral triangle with two other vertices p[a] and p[b];
	// sgn determines the halfplane where the third vertex lies
	point getEquilateralPoint(int a, int b, int sgn) {
		return p[a] + (p[b] - p[a]).rotate(cos60, sgn * sin60);
	}

	// the merge phase of the Melzak algorithm
	//
	// returns the length of the full local Steiner tree on points p[0], ..., p[n - 1] for the given topology,
	// or infinity, if the such tree does not exist;
	//
	// if n = 3, the only possible topology is the following: Steiner point 3 is connected with points 0, 1, 2;
	//
	// if n = 4, the topology is defined by arguments v, a, b:
	// Steiner point 4 is connected to points 0, v, 5; Steiner point 5 is connected to points a, b, 4
	interval merge(int n, int v = -1, int a = -1, int b = -1) {
		if (n == 3) {
			line l(p[1], p[2]);
			
			interval res = inf;
			
			for (int sgn = -1; sgn <= 1; sgn += 2) {
				p[3] = getEquilateralPoint(1, 2, sgn);
				// the equilateral point and p[0] should lie on different sides of l
				if (l.getSide(p[0]) * l.getSide(p[3]) == 1) {
					continue;
				}
				interval length = (p[0] - p[3]).len();
				if (Inf(length) < Inf(res) && reconstruct(3, 1, 2, 0)) {
					res = mymin(res, length);
				}
			}
			return res;
		} else {
			assert(n == 4);

			interval res = inf;

			for (int sgn = -1; sgn <= 1; sgn += 2) {
				for (int sgn2 = -1; sgn2 <= 1; sgn2 += 2) {
					p[5] = getEquilateralPoint(a, b, sgn);
					p[4] = getEquilateralPoint(v, 5, sgn2);

					line l(p[v], p[5]);

					// the equilateral point and p[0] should lie on different sides of l
					if (l.getSide(p[0]) * l.getSide(p[4]) == 1) {
						continue;
					}
					interval length = (p[0] - p[4]).len();
					if (Inf(length) < Inf(res) && reconstruct(4, v, 5, 0) && reconstruct(5, a, b, 4)) {
						res = mymin(res, length);
					}
				}
			}

			return res;
		}
	}

	// find the length of Steiner tree on points inp[0], inp[1], inp[2], inp[3]
	interval solve(point inp[]) {
		// minLength[i] is the length of Steiner tree for the set of points {inp[0], inp[1], inp[2], inp[3]} \ {inp[i]}
		interval minLength[4];

		interval res;

		for (int i = 0; i < 4; ++i) {
			for (int j = 0, z = 0; j < 4; ++j) {
				if (i == j) {
					continue;
				}
				p[z++] = inp[j];
			}
			// {p[0], p[1], p[2]} = {inp[0], inp[1], inp[2], inp[3]} \ {inp[i]}

			// the length of Steiner tree on points p[0], p[1], p[2] is either the length of full local Steiner tree on these points,
			minLength[i] = merge(3);

			// or the sum of two distances between these points
			minLength[i] = mymin(minLength[i], (p[0] - p[1]).len() + (p[0] - p[2]).len());
			minLength[i] = mymin(minLength[i], (p[0] - p[1]).len() + (p[1] - p[2]).len());
			minLength[i] = mymin(minLength[i], (p[0] - p[2]).len() + (p[1] - p[2]).len());
		}

		for (int i = 0; i < 4; ++i) {
			p[i] = inp[i];
		}

		// the length of Steiner tree on all four points
		// is either the length of full local Steiner tree on these points (there are three possible topologies),
		res = mymin(merge(4, 1, 2, 3), mymin(merge(4, 2, 1, 3), merge(4, 3, 1, 2)));

		// or minLength[i] plus the distance from p[i] to the closest point (for some i)
		res = mymin(res, minLength[0] + mymin((p[0] - p[1]).len(), mymin((p[0] - p[2]).len(), (p[0] - p[3]).len())));
		res = mymin(res, minLength[1] + mymin((p[1] - p[0]).len(), mymin((p[1] - p[2]).len(), (p[1] - p[3]).len())));
		res = mymin(res, minLength[2] + mymin((p[2] - p[0]).len(), mymin((p[2] - p[1]).len(), (p[2] - p[3]).len())));
		res = mymin(res, minLength[3] + mymin((p[3] - p[0]).len(), mymin((p[3] - p[1]).len(), (p[3] - p[2]).len())));
		return res;
	}
};

#endif /* MELZAK_H */
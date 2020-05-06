#ifndef MELZAK_H
#define MELZAK_H

#include "interval.hpp"
#include "imath.hpp"
#include "geometry.h"
#include <math.h>
#include <assert.h>

using namespace std;
using namespace cxsc;

// Description of the Melzak algorithm can be found, for example, in "The Steiner Tree Problem" by F.K. Hwang, D.S. Richards, P. Winter.
// In particular, the following terminology from this book is used in the comments:
// -- the Steiner tree;
// -- the Steiner point;
// -- the Steiner minimal tree (SMT);
// -- the full Steiner tree (FST);
// -- the equilateral point.

// This is the implementation of the Melzak algorithm with additional assumption that the number of given points is 4.

// Note that the reconstruction phase is considered successful if the result of some checks could not be determined
// (that usually happens if some of the points are very close to each other).
// It means that the algorithm can return the length smaller that the real length of the SMT on given points.
// This is fine since we are interested only on the lower bound on this length.

// FST contains 6 points (4 given points + 2 Steiner points)
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

		// find the intersection of the circle and the line which goes through p[s] and p[v]
		point res[2];
		intersect(line(p[s], p[v]), o, r, res);

		// l is the line which goes through p[a] and p[b]
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

	// find the equilateral point e, such that e and p[v] are on the opposite sides of the line going through p[a] and p[b];
	// store the found point in p[s]
	void findEquilateralPoint(int s, int a, int b, int v) {
		// l is the line which goes through p[a] and p[b]
		line l(p[a], p[b]);

		// check the equilateral point on one side of l
		point e = p[a] + (p[b] - p[a]).rotate(cos60, sin60);

		if (l.getSide(e) == l.getSide(p[v])) {
			// if e and p[v] are on the same side of l, check the equilateral point on the other side of l
			e = p[a] + (p[b] - p[a]).rotate(cos60, -sin60);
		}

		// check that e and p[v] are on the opposite sides of l
		assert(l.getSide(e) * l.getSide(p[v]) == -1);

		p[s] = e;
	}

	// the merge phase of the Melzak algorithm
	//
	// returns the length of the FST on points p[0], ..., p[n - 1] for the given topology,
	// or infinity, if the FST does not exist;
	//
	// if n = 3, the only possible topology is the following: Steiner point 3 is connected with 0, 1, 2;
	//
	// if n = 4, the topology is defined by arguments v, a, b:
	// Steiner point 4 is connected to 0, v, 5; Steiner point 5 is connected to a, b, 4
	interval merge(int n, int v = -1, int a = -1, int b = -1) {
		if (n == 3) {
			findEquilateralPoint(3, 1, 2, 0);
			interval length = (p[0] - p[3]).len(); // the length of the FST (if it exists)

			return (reconstruct(3, 1, 2, 0) ? length : inf);
		} else {
			assert(n == 4);

			// l is the line which goes through p[a] and p[b]
			line l(p[a], p[b]);

			interval res = inf;

			// check equilateral points on both sides of line l
			for (int sgn = -1; sgn <= 1; sgn += 2) {
				p[5] = p[a] + (p[b] - p[a]).rotate(cos60, sgn * sin60);

				findEquilateralPoint(4, v, 5, 0);
				interval length = (p[0] - p[4]).len(); // the length of the FST (if it exists)

				if (Inf(length) < Inf(res) && reconstruct(4, v, 5, 0) && reconstruct(5, a, b, 4)) {
					res = mymin(res, length);
				}
			}

			return res;
		}
	}

	// find the length of the SMT on points inp[0], inp[1], inp[2], inp[3]
	interval solve(point inp[]) {
		// minLength[i] is the length of the SMT for the set of points {inp[0], inp[1], inp[2], inp[3]} \ {inp[i]}
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

			// the length of the SMT on points p[0], p[1], p[2] is either the length of the FST on these points,
			minLength[i] = merge(3);

			// or the sum of two distances between these points
			minLength[i] = mymin(minLength[i], (p[0] - p[1]).len() + (p[0] - p[2]).len());
			minLength[i] = mymin(minLength[i], (p[0] - p[1]).len() + (p[1] - p[2]).len());
			minLength[i] = mymin(minLength[i], (p[0] - p[2]).len() + (p[1] - p[2]).len());
		}

		for (int i = 0; i < 4; ++i) {
			p[i] = inp[i];
		}

		// the length of the SMT on all four points
		// is either the length of the FST on these points (there are three possible topologies),
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
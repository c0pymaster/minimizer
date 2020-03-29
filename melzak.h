#ifndef MELZAK_H
#define MELZAK_H

#include "geometry.h"
#include <math.h>
#include <assert.h>

#define eprintf(...) fprintf(stderr, __VA_ARGS__), fflush(stderr)
#define sz(a) ((int) (a).size())

using namespace std;

// Description of the Melzak algorithm can be found, for example, in "The Steiner Tree Problem" by F.K. Hwang, D.S. Richards, P. Winter.
// In particular, the following terminology from this book is used in the comments:
// -- the Steiner tree;
// -- the Steiner point;
// -- the Steiner minimal tree (SMT);
// -- the full Steiner tree (FST);
// -- the equilateral point.

// This is the implementation of the Melzak algorithm with additional assumption that the number of points is 4.

// FST contains 6 points (4 initial points + 2 Steiner points)
const int maxn = 6;

struct Melzak {
	point p[maxn];

	// the reconstruction step of the Melzak algorithm;
	// returns true if the reconstruction is successfull
	bool reconstruct(int s, int a, int b, int v) {
		// construct the circle through points p[a], p[b], p[s]
		point o;
		long double r;
		getCircle(p[a], p[b], p[s], o, r);

		// find the intersection of the circle and the line which goes through p[s] and p[v]
		point res[2];
		intersect(line(p[s], p[v]), o, r, res);

		// check that the intersection is not empty (it should be nonempty since p[s] is both on the line and on the circle)
		assert(intersect(line(p[s], p[v]), o, r, res));

		// if one of the intersection points is on the inside of both the segment between p[s] and p[v], and the smaller arc between p[a] and p[b],
		// then this point is the desired Steiner point; otherwise, the reconstruction phase failed
		bool found = false;
		for (int i = 0; i < 2; ++i) {
			if (!(p[s] == res[i]) && insideSeg(p[s], p[v], res[i]) && insideArc(p[a] - o, p[b] - o, res[i] - o)) {
				found = true;

				// check that the length of the tree does not change during the reconstruction step
				assert(abs((p[s] - res[i]).len() - (p[a] - res[i]).len() - (p[b] - res[i]).len()) < eps);

				p[s] = res[i];
				break;
			}
		}

		return found;
	}

	// try to find the equilateral point e for points p[a], p[b],
	// such that e and p[v] are on the opposite sides of the line going through p[a] and p[b];
	// if such e exists, store it in p[s]
	bool findEquilateralPoint(int s, int a, int b, int v) {
		if (p[a] == p[b]) {
			// if p[a] = p[b], e does not exist
			return false;
		}

		// l is the line which goes through p[a] and p[b]
		line l(p[a], p[b]);

		if (l.getSide(p[v]) == 0) {
			// if p[v] is on l, e does not exist
			return false;
		}

		// checking the point on one side of l
		point e = p[a] + (p[b] - p[a]).rotate(cos60, sin60);
		assert(equilateralTriangle(p[a], p[b], e));

		if (l.getSide(e) == l.getSide(p[v])) {
			// if e and p[v] are on the same side of l, checking the point on the other side of l
			e = p[a] + (p[b] - p[a]).rotate(cos60, -sin60);
			assert(equilateralTriangle(p[a], p[b], e));
		}

		// check that e and p[v] are on the opposite sides of l
		assert(l.getSide(e) * l.getSide(p[v]) == -1);

		p[s] = e;

		return true;
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
	long double merge(int n, int v = -1, int a = -1, int b = -1) {
		if (n == 3) {
			// if the equilateral point exists, call the reconstruction step;
			// otherwise there is no FST for this set of points
			if (findEquilateralPoint(3, 1, 2, 0)) {
				// the length of the FST (if it exists)
				long double length = (p[0] - p[3]).len();
				
				if (reconstruct(3, 1, 2, 0)) {
					// the reconstruction is successfull
					return length;
				}
			}

			return inf;
		} else {
			assert(n == 4);

			// l is the line which goes through p[a] and p[b]
			line l(p[a], p[b]);

			long double res = inf;

			// checking equilateral points on both sides of line l
			for (int sgn = -1; sgn <= 1; sgn += 2) {
				p[5] = p[a] + (p[b] - p[a]).rotate(cos60, sgn * sin60);
				assert(equilateralTriangle(p[a], p[b], p[5]));

				// if the second equilateral point exists, call reconstruction steps;
				// otherwise there is no FST for this topology and this choice of the first equilateral point
				if (findEquilateralPoint(4, v, 5, 0)) {
					// the length of the FST (if it exists)
					long double length = (p[0] - p[4]).len();

					if (length < res && reconstruct(4, v, 5, 0) && reconstruct(5, a, b, 4)) {
						// the reconstruction is successfull
						res = length;
					}
				}
			}

			return res;
		}
	}

	// find the length of the SMT on points inp[0], inp[1], inp[2], inp[3]
	long double solve(point inp[]) {
		// minLength[i] is the length of the SMT for the set of points {inp[0], inp[1], inp[2], inp[3]} \ {inp[i]}
		long double minLength[maxn];

		for (int i = 0; i < 4; ++i) {
			for (int j = 0, z = 0; j < 4; ++j) {
				if (i == j) {
					continue;
				}
				p[z++] = inp[j];
			}
			// {inp[0], inp[1], inp[2], inp[3]} \ {inp[i]} = {p[0], p[1], p[2]}

			// the length of the SMT on points p[0], p[1], p[2] is either the length of the FST on these points,
			minLength[i] = merge(3);

			// or the sum of two distances between these points
			minLength[i] = min(minLength[i], (p[0] - p[1]).len() + (p[0] - p[2]).len());
			minLength[i] = min(minLength[i], (p[0] - p[1]).len() + (p[1] - p[2]).len());
			minLength[i] = min(minLength[i], (p[0] - p[2]).len() + (p[1] - p[2]).len());
		}

		for (int i = 0; i < 4; ++i) {
			p[i] = inp[i];
		}

		// the length of the SMT on all four points
		// is either the length of the FST on these points (there are three possible topologies),
		long double res = min(merge(4, 1, 2, 3), min(merge(4, 2, 1, 3), merge(4, 3, 1, 2)));

		// or minLength[i] plus the distance from p[i] to the closest point (for some i)
		res = min(res, minLength[0] + min((p[0] - p[1]).len(), min((p[0] - p[2]).len(), (p[0] - p[3]).len())));
		res = min(res, minLength[1] + min((p[1] - p[0]).len(), min((p[1] - p[2]).len(), (p[1] - p[3]).len())));
		res = min(res, minLength[2] + min((p[2] - p[0]).len(), min((p[2] - p[1]).len(), (p[2] - p[3]).len())));
		res = min(res, minLength[3] + min((p[3] - p[0]).len(), min((p[3] - p[1]).len(), (p[3] - p[2]).len())));

		return res;
	}
};

#endif /* MELZAK_H */
#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "interval.hpp"
#include "imath.hpp"
#include <math.h>
#include <assert.h>
#include <algorithm>
#include <vector>

using namespace std;
using namespace cxsc;

const interval pi = Pi();
const interval sin60 = sqrt(interval(3)) / interval(2);
const interval cos60 = interval(1) / interval(2);
const interval inf = interval(1e18);

interval mymin(const interval & a, const interval & b) {
	return interval(min(Inf(a), Inf(b)), min(Sup(a), Sup(b)));
}

interval mymax(const interval & a, const interval & b) {
	return interval(max(Inf(a), Inf(b)), max(Sup(a), Sup(b)));
}

// return the sign of the real number x;
// if the sign cannot be determined, return 0
int getSign(const interval & x) {
	if (Inf(x) <= 0 and 0 <= Sup(x)) { // sign cannot be determined
		return 0;
	} else if (0 < Inf(x)) {
		return 1;
	} else {
		assert(Sup(x) < 0);
		return -1;
	}
}

// a point or a vector on a plane
struct point {
	interval x, y;

	point(): x(), y() {}

	point(const interval & x, const interval & y): x(x), y(y) {}

	// return the difference of two vectors
	point operator -(const point & p) const {
		return point(x - p.x, y - p.y);
	}

	// return the sum of two vectors
	point operator +(const point & p) const {
		return point(x + p.x, y + p.y);
	}

	// return the dot product of two vectors
	interval operator *(const point & p) const {
		return x * p.x + y * p.y;
	}

	// return the skew product of two vectors
	interval operator ^(const point & p) const {
		return x * p.y - y * p.x;
	}

	// return the vector, multiplied by the constant
	point operator *(const interval & m) const {
		return point(x * m, y * m);
	}

	// return the squared length of the vector
	interval slen() const {
		return power(x, 2) + power(y, 2);
	}

	// return the length of the vector
	interval len() const {
		return sqrt(slen());
	}

	// return the normalized vector
	point normalize() const {
		interval l = len();

		// check that the length of the vector is not zero
		assert(Inf(l) > 0);

		return point(x / l, y / l);
	}

	// return the vector, rotated by pi / 2 counter-clockwise
	point rotate90() const {
		return point(-y, x);
	}

	// return the vector, rotated by angle radians counter-clockwise, where cos(angle) = cosa, sin(angle) = sina
	point rotate(const interval & cosa, const interval & sina) const {
		return point(x * cosa - y * sina, x * sina + y * cosa);
	}

	// return the vector, rotated by angle radians counter-clockwise
	point rotate(const interval & angle) const {
		return rotate(cos(angle), sin(angle));
	}
};

// debug output
void eprint(const vector<point> & p) {
	for (int i = 0; i < (int) p.size(); ++i) {
		cerr << "(" << p[i].x << ", " << p[i].y << ")\n";
	}
	cerr.flush();
}

// a line on a plane, described by a vector equation o + v * t, where t is the parameter
struct line {
	point o, v;

	line(): o(), v() {}

	// construct the line passing through two points a and b
	line(const point & a, const point & b) {
		o = a;
		v = (b - a).normalize();
	}

	// return the signed distance from point p to the line (positive in one half-plane, negative in the other)
	interval dist(const point & p) const {
		return (p - o) ^ v;
	}

	// return the half-plane identifier (+1 or -1);
	// return 0, if the half-plane cannot be determined
	int getSide(const point & p) const {
		return getSign(dist(p));
	}
};

// intersect line l and the circle with center o and radius r;
// assume the intersection is not empty
void intersect(const line & l, const point & o, const interval & r, point res[2]) {
	point n = l.v.rotate90();
	// n is the normal vector of line l
	
	interval d = l.dist(o);
	interval x = sqrt(mymax(interval(0), power(r, 2) - power(d, 2)));
	for (int it = 0; it < 2; ++it) {
		res[it] = o + n * d + l.v * (it ? -x : x);

		// the found points should lie both on the circle and on the line
		assert(getSign(r - (o - res[it]).len()) == 0 && getSign(l.dist(res[it])) == 0);
	}
}

// intersect lines a and b, assume they do not coincide and are not parallel
point intersect(const line & a, const line & b) {
	interval coeff = (a.v ^ b.v);

	// check that the lines do not coincide and are not parallel
	assert(getSign(coeff) != 0);

	interval ka = (b.o ^ b.v) / coeff;
	interval kb = (a.v ^ a.o) / coeff;
	point res = a.v * ka + b.v * kb;

	// the found point should lie on both given lines
	assert(getSign(a.dist(res)) == 0 && getSign(b.dist(res)) == 0);

	return res;
}

// get the perpendicular bisector to the segment ab
line getPerpendicularBisector(const point & a, const point & b) {
	line l;
	l.o = (a + b) * interval(0.5);
	l.v = (b - a).rotate90().normalize();
	return l;
}

// get the circle that passes through points a, b, c (assume that the circle exists)
void getCircle(const point & a, const point & b, const point & c, point & o, interval & r) {
	o = intersect(getPerpendicularBisector(a, b), getPerpendicularBisector(a, c));
	r = (a - o).len();

	// the found circle should contain all three given points
	assert(getSign(r - (a - o).len()) == 0 && getSign(r - (b - o).len()) == 0 && getSign(r - (c - o).len()) == 0);
}

#endif /* GEOMETRY_H */
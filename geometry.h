#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <math.h>
#include <assert.h>
#include <algorithm>
#include <vector>

#define eprintf(...) fprintf(stderr, __VA_ARGS__), fflush(stderr)
#define sz(a) ((int) (a).size())

using namespace std;

const long double inf = 1e18;
const long double eps = 1e-8;

const long double pi = M_PI;
const long double sin60 = sqrt((long double) 3.0) / 2;
const long double cos60 = 0.5;

// get sign of the real number x
int getSign(const long double & x) {
	if (abs(x) < eps) {
		return 0;
	}
	return ((x > 0) ? 1 : -1);
}

// a point or a vector on a plane
struct point {
	long double x, y;

	point(): x(), y() {}

	point(const long double & x, const long double & y): x(x), y(y) {}

	// return the difference of two vectors
	point operator -(const point & p) const {
		return point(x - p.x, y - p.y);
	}

	// return the sum of two vectors
	point operator +(const point & p) const {
		return point(x + p.x, y + p.y);
	}

	// return the dot product of two vectors
	long double operator *(const point & p) const {
		return x * p.x + y * p.y;
	}

	// return the skew product of two vectors
	long double operator ^(const point & p) const {
		return x * p.y - y * p.x;
	}

	// return the vector, multiplied by the constant
	point operator *(const long double & m) const {
		return point(x * m, y * m);
	}

	// return the squared length of the vector
	long double slen() const {
		return x * x + y * y;
	}

	// return the length of the vector
	long double len() const {
		return sqrt(slen());
	}

	// return the normalized vector
	point normalize() const {
		long double l = len();

		// check that the length of the vector is not zero
		assert(l > eps);

		return point(x / l, y / l);
	}

	// return the vector, rotated by pi / 2 counter-clockwise
	point rotate90() const {
		return point(-y, x);
	}

	// return the vector, rotated by angle radians counter-clockwise, where cos(angle) = cosa, sin(angle) = sina
	point rotate(const long double & cosa, const long double & sina) const {
		return point(x * cosa - y * sina, x * sina + y * cosa);
	}

	// return the vector, rotated by angle radians counter-clockwise
	point rotate(const long double & angle) const {
		return rotate(cos(angle), sin(angle));
	}

	// compare two points by their x-coordinates, and then by their y-coordinates
	bool operator <(const point & p) const {
		if (abs(x - p.x) > eps) {
			return x < p.x;
		}
		if (abs(y - p.y) > eps) {
			return y < p.y;
		}
		return false;
	}

	// check if two points coincide
	bool operator ==(const point & p) const {
		return abs(x - p.x) < eps && abs(y - p.y) < eps;
	}
};

// debug output
void eprint(const vector<point> & p) {
	for (int i = 0; i < sz(p); ++i) {
		eprintf("(%.3Lf,%.3Lf)%c", p[i].x, p[i].y, " \n"[i == sz(p) - 1]);
	}
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
	long double dist(const point & p) const {
		return (p - o) ^ v;
	}

	// return 0 if point p is on the line; otherwise return the half-plane identifier (+1 or -1);
	// this function is used to check if two points are in the same half-plane
	int getSide(const point & p) const {
		return getSign(dist(p));
	}
};

// points a, b, c are on the same line, a != b; check if c is the interior of the segment ab
bool insideSeg(const point & a, const point & b, const point & c) {
	// check that a, b, c are on the same line
	assert(abs((b - a) ^ (c - a)) < eps);

	if (a < b) {
		return a < c && c < b;
	} else {
		// since a != b, either a < b or b < a
		assert(b < a);

		return b < c && c < a;
	}
}

// a, b, c are points on the circle with the center at (0, 0)
// check that segment ab is not a diameter, and that point c is inside the interior of the smaller arc between points a and b
//
// used in the reconstruction phase of the Melzak algorithm;
// the reconstruction is successful only if the angle between a and b is 2 * pi / 3 (so ab is not a diameter)
bool insideArc(const point & a, const point & b, const point & c) {
	// check that a, b, c have the same length
	assert(abs(a.slen() - b.slen()) < eps && abs(a.slen() - c.slen()) < eps);

	int sab = getSign(a ^ b), sac = getSign(a ^ c), scb = getSign(c ^ b);
	return (sab != 0 && sab == sac && sab == scb);
}

// intersect line l and the circle with center o and radius r
bool intersect(const line & l, const point & o, const long double & r, point res[2]) {
	point n = l.v.rotate90();
	// n is the normal vector of line l
	
	long double d = l.dist(o);
	if (abs(d) > r + eps) {
		// the circle and the line do not intersect
		return false;
	}
	long double x = sqrt(max((long double) 0, r * r - d * d));
	for (int it = 0; it < 2; ++it) {
		res[it] = o + n * d + l.v * (it ? -x : x);
	}
	return true;
}

// intersect lines a and b, assume they do not coincide and are not parallel
point intersect(const line & a, const line & b) {
	long double coeff = (a.v ^ b.v);

	// check that the lines do not coincide and are not parallel
	assert(abs(coeff) > eps);

	long double ka = (b.o ^ b.v) / coeff;
	long double kb = (a.v ^ a.o) / coeff;
	point res = a.v * ka + b.v * kb;

	// check that the intersection point is lying on both lines
	assert(abs(a.dist(res)) < eps && abs(b.dist(res)) < eps);

	return res;
}

// get the perpendicular bisector to the segment ab
line getPerpendicularBisector(const point & a, const point & b) {
	// check that a != b
	assert(!(a == b));

	line l;
	l.o = (a + b) * 0.5;
	l.v = (b - a).rotate90().normalize();
	return l;
}

// get the circle that passes through points a, b, c (assume that the circle exists)
void getCircle(const point & a, const point & b, const point & c, point & o, long double & r) {
	// check that a, b, c are pairwise distinct
	assert(!(a == b || a == c || b == c));

	o = intersect(getPerpendicularBisector(a, b), getPerpendicularBisector(a, c));
	r = (a - o).len();

	// check that |oa| = |ob| = |oc|
	assert(abs(r - (b - o).len()) < eps && abs(r - (c - o).len()) < eps);
}

//check if the triangle on points a, b, c is equilateral
bool equilateralTriangle(const point & a, const point & b, const point & c) {
	return (abs((a - b).slen() - (a - c).slen()) < eps) && (abs((a - b).slen() - (b - c).slen()) < eps);
}


#endif /* GEOMETRY_H */
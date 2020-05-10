/* Include file for weil.cc
 * Defines the Point class, and the Fp2 class
 * Ben Lynn
 */
/*
Copyright (C) 2001 Benjamin Lynn (blynn@cs.stanford.edu)

This file is part of the Stanford short signature system.

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/
#ifndef WEIL_H
#define WEIL_H

#include "f3kl.h"

typedef F3kl Fp2;

//a hack right now; starts with z = 0, can't actually set z = 1,
//use put_O to set to O = (0, 1, 0)
//allows all 2D points and one point at infinity
class Point {
public:
	Fp2 x, y;
	bool infinity;
	void put(Fp2 a, Fp2 b) { x=a; y=b; infinity = false; }
	void put_O(void) { set(x); set(y); infinity = true; }
	void set_finite(void) { infinity = false; }
	void random(void);
	bool is_O(void) { return (infinity); } // && IsZero(x) && IsOne(y)); }
	Point () { infinity = false; }

	friend ostream& operator<<(ostream& s, const Point& a) {
		if (a.infinity) {
			s << "O";
		} else {
			s << a.x << " " << a.y;
		}
		return s;
	}

/*
	friend istream& operator>>(istream& s, Point& a) {
		if (s.peek() == (int) 'O') {
			a.put_O();
			return s;
		}
		a.set_finite();
		s >> a.x;
		s >> a.y;
		return s;
	}
	*/
};

long operator !=(const Point& P, const Point& Q);
long operator ==(const Point& P, const Point& Q);

void random_point(Point& P);
void group_plus(Point &R, const Point &P, const Point &Q);
void group_times(Point &R, const ZZ &a, const Point &P);
Fp2 tate_pairing(const Point& P, const Point& Q);
Fp2 weil_pairing(const Point& P, const Point& Q);
Fp2 tate_power(Fp2& input);
Fp2 ehat(const Point& P, const Point& Q);
void get_line(const Point& nom, const Point& denom, const Point& P, const Point& Q);
void get_vertical(const Point& nom, const Point& denom, const Point& P);
//Fp2 eval(const Fp2& a, const Fp2& b, const Fp2& c, const Point& P);
void weil_init_order(int l, const ZZ&q, const ZZ&multiplier);
void weil_init_nqr(int);
void weil_init_curve(Fp2, Fp2, Fp2);
bool y_from_x(Fp2& y, const Fp2& x);
void force_subgroup(Point& P);

#endif WEIL_H

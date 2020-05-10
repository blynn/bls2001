/* f3k.h
  defines operations for elements in F_3^k
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
#ifndef F3K_H
#define F3K_H

#include <NTL/ZZ.h>

extern int discretelog[], antilog[], addtable[];
//from f3ktest.cc
enum { phirminusexponent = 476,
	phirplusexponent = 112,
	phiuexponent = 182 };

enum { f3korder = 728 };

class F3k {
public:
	int b;

	F3k() {
		b = 0;
	}

	inline void put(int a) {
		b = discretelog[a];
	}
	F3k& operator=(const F3k& a) {
		b = a.b;
		return *this;
	}

	friend ostream& operator<<(ostream& s, const F3k& a) {
		//s << "[" << a.b1 << " " << a.b2 << "]";
		if (a.b < 0) {
			s << 0;
		} else {
			s << antilog[a.b];
		}
		return s;
	}
};

inline void sqr(F3k&x, const F3k &a) {
	if (a.b < 0) {
		x.b = -(f3korder + 1);
	} else {
		x.b = a.b * 2;
		if (x.b >= f3korder) x.b -= f3korder;
	}
}

inline void mul(F3k&x, const F3k &a, const F3k& b) {
	x.b = a.b + b.b;
	if (x.b >= f3korder) x.b -= f3korder;
	else if (x.b < 0) {
		x.b = -(f3korder + 1);
	}
}

inline void inv(F3k& res, const F3k& a) {
	if (a.b < 0) {
		cerr << "f3k: can't invert 0\n";
		exit(1);
	} else if (a.b > 0) {
		res.b = f3korder - a.b;
	}
}
inline void add(F3k& x, const F3k& a, const F3k& b) {
	//use a+b = a(1+a^-1 b) trick
	if (a.b < 0) {
		x.b = b.b;
		return;
	}
	if (b.b < 0 ) {
		x.b = a.b;
		return;
	}
	int a1;
	a1 = f3korder - a.b;
	int ab;
	ab = a1 + b.b;
	if (ab >= f3korder) ab -= f3korder;
	a1 = addtable[ab];
	x.b = a1 + a.b;
	if (x.b >= f3korder) x.b -= f3korder;
	else if (x.b < 0) {
		x.b = -(f3korder + 1);
	}
}

inline void negate(F3k& res, const F3k& a) {
	//res.b = negatetable[a.b];
	if (a.b < 0) {
		res.b = -(f3korder + 1);
		return;
	}
	res.b = a.b + (f3korder / 2);
	if (res.b >= f3korder) res.b -= f3korder;
}

inline void div(F3k& res, const F3k& a, const F3k& b)
{
	F3k temp;
	inv(temp, b);
	mul(res, a, temp);
}

inline void clear(F3k& x) {
	x.b = -(f3korder + 1);
}
inline void set(F3k& x) {
	x.b = 0;
}

inline bool IsZero(const F3k& a) {
	return (a.b < 0);
}

inline bool IsOne(const F3k& a) {
	return (a.b == 0);
}

inline long operator!=(const F3k& a, const F3k& b) {
	return (a.b != b.b);
}
inline long operator==(const F3k& a, const F3k& b) {
	return (a.b == b.b);
}
#endif F3K_H

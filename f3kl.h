/* f3kl.h
  defines operations for elements in F_3^kl
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
#ifndef F3KL_H
#define F3KL_H

#include "f3k.h"
//#include <string.h> //for memcpy - seems slower!

class F3kl {
public:
	//F3k *poly;
	F3k poly[170];
	int degree;
	int space;

	F3kl() {
		degree = -1;
		space = 170;
		//int i;
		//for (i=0; i<space; i++) {
			//clear(poly[i]);
		//}
		//poly = new F3k[space+1];
	}

	~F3kl() {
		//delete[] poly;
	}

	inline void graduate(void)
	//compute new degree
	{
		int i;
		for (i=degree; i>=0; i--) {
			if (!IsZero(poly[i])) break;
		}
		degree = i;
	}

	inline void SetCoeff(int d, const F3k &a)
	//warning: if a = 0 and d = degree, call
	//graduate() to get new degree
	{
		if (d < space) {
			poly[d] = a;
		} else {
			cerr << "big polys not supported yet!\n";
		}
		if (d > degree) {
			int i;
			for (i=degree+1; i<d; i++) {
				clear(poly[i]);
			}
			degree = d;
		}
	}

	inline void SetCoeff(int d, int a) {
		if (d < space) {
			switch(a) {
				case 0:
					clear(poly[d]);
					break;
				case 1:
					set(poly[d]);
					break;
				case 2:
					set(poly[d]);
					negate(poly[d], poly[d]);
					break;
			}
		} else {
			cerr << "big polys not supported yet!\n";
		}
		if (d > degree) {
			int i;
			for (i=degree+1; i<d; i++) {
				clear(poly[i]);
			}
			degree = d;
		}
	}

	F3kl& operator=(const F3kl& a) {
		if (a.degree > space) {
			cerr << "big polys not supported yet!\n";
		} else {
			int i;
			degree = a.degree;
			for (i=0; i<=degree; i++) {
				poly[i] = a.poly[i];
			}
		}
		return *this;
	}

	friend ostream& operator<<(ostream& s, const F3kl& a) {
		s << "[";
		int i;
		for (i=0; i<=a.degree; i++) {
			s << a.poly[i];
			if (i < a.degree) s << " ";
		}
		s << "]";
		return s;
	}

};

void f3kl_init(int l);
extern F3kl minpoly;
extern int minpolydegree;
extern F3kl f3klscratch0;
extern F3kl f3klscratch1;

inline void set(F3kl& x) {
	x.degree = 0;
	set(x.poly[0]);
}

inline void clear(F3kl& x) {
	x.degree = -1;
}

inline bool IsOne(const F3kl& x) {
	return x.degree == 0 && IsOne(x.poly[0]);
}
inline bool IsZero(const F3kl& x) {
	return x.degree == -1;
}

inline void add(F3kl& x, const F3kl& a, const F3kl& b) {
	int i, d;
	/*XXX: check space
		cerr << "can't handle this!\n";
		exit(1);
	*/
	if (a.degree > b.degree) {
		d = a.degree;
		for (i=0; i<=b.degree; i++) {
			add(x.poly[i], a.poly[i], b.poly[i]);
		}
		if (&x != &a) {
			/*
			memcpy(&(x.poly[i]), &(a.poly[i]), sizeof(F3k) * (a.degree - i + 1));
			*/
			for (; i<=a.degree; i++) {
				x.poly[i] = a.poly[i];
			}
			/* pointer arithmetic; doesn't seem to speed up much
			F3k *p, *pmax;
			const F3k *q;
			pmax = &x.poly[a.degree];
			p = &x.poly[i];
			q = &a.poly[i];
			for (;p<=pmax; p++, q++) {
				*p = *q;
			}
			*/
		}

	} else {
		d = b.degree;
		for (i=0; i<=a.degree; i++) {
			add(x.poly[i], a.poly[i], b.poly[i]);
		}
		if (&x != &b) {
			/*
			memcpy(&(x.poly[i]), &(b.poly[i]), sizeof(F3k) * (b.degree - i + 1));
			*/
			for (; i<=b.degree; i++) {
				x.poly[i] = b.poly[i];
			}
		}
	}
	for (i=d; i>=0; i--) {
		if (!IsZero(x.poly[i])) break;
	}
	x.degree = i;
}

void addx(F3kl& x, const F3kl& a, const F3kl& b, int pwr);

inline void negate(F3kl& res, const F3kl& a) {
	if (res.space < a.degree) {
		cerr << "can't handle this yet!\n";
		exit(1);
	}
	res.degree = a.degree;
	int i;
	for (i=0; i<=res.degree; i++) {
		negate(res.poly[i], a.poly[i]);
	}
}

inline void sub(F3kl& x, const F3kl& a, const F3kl& b) {
	negate(f3klscratch0, b);
	add(x, f3klscratch0, a);
}

inline void subx(F3kl& x, const F3kl& a, const F3kl& b, int pwr) {
	negate(f3klscratch0, b);
	addx(x, a, f3klscratch0, pwr);
}

void mul(F3kl&x, const F3k &a, const F3kl& b);
void mul(F3kl&x, const F3kl &a, const F3kl& b);

inline void sqr(F3kl&x, const F3kl& a) {
	//XXX: are there optimizations?
	mul(x, a, a);
}

//internal use only
inline void realmul(F3kl&x, const F3kl &a, const F3kl& b)
//multiplication without reducing
{
	int i;
	F3kl& w1 = f3klscratch0;
	F3kl& w2 = f3klscratch1;
	w1 = a;
	clear(x);
	for (i=0; i<=b.degree; i++) {
		mul(w2, b.poly[i], w1);
		add(x, x, w2);
		int j;
		//mult by x
		w1.degree++;
		for (j=w1.degree; j>0; j--) {
			w1.poly[j] = w1.poly[j - 1];
		}
		clear(w1.poly[0]);
	}
}

void divrem(F3kl& q, F3kl& r, const F3kl& a, const F3kl& b);

void egcd(F3k &d, F3kl &x, F3kl &y, const F3kl& a, const F3kl &b);

inline void inv(F3kl& res, const F3kl& a) {
	F3k d;
	
	egcd(d, f3klscratch0, res, minpoly, a);
	inv(d, d);
	mul(res, d, res);
}

inline void div(F3kl& q, const F3kl& a, const F3kl& b)
{
	inv(f3klscratch0, b);
	mul(q, a, f3klscratch0);
}

inline long operator!=(const F3kl& a, const F3kl& b) {
	if (a.degree != b.degree) return true;
	int i;
	for (i=0; i<a.degree; i++) {
		if (a.poly[i] != b.poly[i]) return true;
	}
	return false;
}

inline long operator==(const F3kl& a, const F3kl& b) {
	if (a.degree != b.degree) return false;
	int i;
	for (i=0; i<a.degree; i++) {
		if (a.poly[i] != b.poly[i]) return false;
	}
	return true;
}

void power(F3kl& res, const F3kl& x, const ZZ& n);

void initminpoly(int l);

#endif F3KL_H

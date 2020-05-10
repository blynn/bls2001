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
#include "f3kl.h"

F3kl minpoly;
int minpolydegree;
F3kl replacepoly[729];
F3kl f3klscratch0;
F3kl f3klscratch1;

void initminpoly(int l) {

	switch(l) {
		case 79:
			minpoly.SetCoeff(0, 1);
			minpoly.SetCoeff(3, 2);
			minpoly.SetCoeff(4, 1);
			minpoly.SetCoeff(79, 1);
			break;
		case 97:
			minpoly.SetCoeff(0, 1);
			minpoly.SetCoeff(1, 2);
			minpoly.SetCoeff(2, 2);
			minpoly.SetCoeff(3, 2);
			minpoly.SetCoeff(97, 1);
			break;
		case 121:
			minpoly.SetCoeff(0, 1);
			minpoly.SetCoeff(1, 2);
			minpoly.SetCoeff(121, 1);
			break;
		case 149:
			minpoly.SetCoeff(0, 1);
			minpoly.SetCoeff(1, 2);
			minpoly.SetCoeff(3, 1);
			minpoly.SetCoeff(4, 1);
			minpoly.SetCoeff(6, 1);
			minpoly.SetCoeff(149, 1);
			break;
		case 163:
			minpoly.SetCoeff(0, 1);
			minpoly.SetCoeff(2, 1);
			minpoly.SetCoeff(3, 2);
			minpoly.SetCoeff(4, 1);
			minpoly.SetCoeff(5, 1);
			minpoly.SetCoeff(163, 1);
			break;
		case 167:
			minpoly.SetCoeff(0, 1);
			minpoly.SetCoeff(1, 1);
			minpoly.SetCoeff(2, 2);
			minpoly.SetCoeff(3, 2);
			minpoly.SetCoeff(167, 1);
			break;
	}

	minpolydegree = minpoly.degree;
	replacepoly[0] = minpoly;
	replacepoly[0].SetCoeff(minpolydegree, 0);
	negate(replacepoly[0], replacepoly[0]);
	replacepoly[0].graduate();
	int i;
	for (i=1; i<729; i++) {
		F3k m;
		m.b = i;
		mul(replacepoly[i], m, replacepoly[0]);
	}
}

enum { windowsize = 4,
	windowsizepower = 7		//this is 2^windowsize - 1
};

static int intpower2[windowsize];

void f3kl_init(int l)
{
	int i;
	intpower2[0] = 1;
	for (i=1; i<windowsize; i++) intpower2[i] = intpower2[i-1] * 2;
	initminpoly(l);
}

void power(F3kl&res, const F3kl&x, const ZZ&n)
{
	int i;

	set(res);

	//use sliding-window method
	F3kl g[2*windowsizepower+2];
	g[1] = x;
	mul(g[2], x, x);

	for (i=1; i<=windowsizepower; i++) {
		mul(g[2*i+1], g[2*i-1], g[2]);
	}
	int m = NumBits(n) - 1;

	while(m>=0) {
		if (!bit(n, m)) {
			mul(res, res, res);
			m--;
		} else {
			int l;
			l = m - windowsize + 1;
			if (l < 0) l = 0;
			for (; l<m; l++) {
				if (bit(n, l)) break;
			}
			int j = 1;
			int k;
			mul(res, res, res);
			for (k=l+1; k<=m; k++) {
				if (bit(n, k)) j += intpower2[k-l];
				mul(res, res, res);
			}
			mul(res, res, g[j]);
			m = l-1;
		}
	}
}

void mul(F3kl&x, const F3k &a, const F3kl& b)
{
	if (x.space < b.degree) {
		cerr << "can't handle this yet!\n";
		exit(1);
	}
	x.degree = b.degree;
	int i;
	for (i=0; i<=x.degree; i++) {
		mul(x.poly[i], a, b.poly[i]);
	}
}

void mul(F3kl&res, const F3kl &a, const F3kl& b)
//replace with Horner's rule?
{
	int i;
	F3kl w1, w2, x;
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
		if (w1.degree == minpolydegree) {
			//reduce by minpoly
			//mul(w2, w1.poly[minpolydegree], replacepoly);
			//add(w1, w1, w2);
			int m;
			m = w1.poly[minpolydegree].b;
			clear(w1.poly[minpolydegree]);
			add(w1, replacepoly[m], w1);
		}
	}
	res = x;
}

void divrem(F3kl& q, F3kl& r, const F3kl& a, const F3kl& b)
{
	int da, db;
	da = a.degree;
	db = b.degree;

	if (db < 0) {
		cerr << "f3kl: division by zero\n";
		exit(1);
	}
	r = a;
	clear(q);
	if (da < db) {
		return;
	}

	F3kl res, temp;
	F3k c;
	F3k lb;
	inv(lb, b.poly[db]);
	int i;
	int dr;
	for (;;) {
		dr = r.degree;
		i = dr - db;
		if (dr < db) break;
		mul(c, r.poly[dr], lb);
		res.SetCoeff(i, c);
		mul(temp, c, b);
		subx(r, r, temp, dr - db);
	}
	q = res;
}

void egcd(F3k &d, F3kl &x, F3kl &y, const F3kl& a, const F3kl &b)
{
	if (IsZero(b)) {
		//cout << "bottom: " << a << endl;
		d = a.poly[0];
		set(x);
		clear(y);
		return;
	}
	F3kl x1, y1;
	F3kl q, r;
	divrem(q, r, a, b);
	egcd(d, x1, y1, b, r);
	x = y1;
	mul(y1, y1, q);
	sub(y, x1, y1);
}

void addx(F3kl& x, const F3kl& a, const F3kl& b, int pwr)
//x = a + b*x^pwr
{
	int i, d;
	/*XXX: check space
		cerr << "can't handle this!\n";
		exit(1);
	*/
	int bpwr = b.degree + pwr;

	if (a.degree >= bpwr) {
		d = a.degree;
		for (i=0; i<pwr; i++) {
			x.poly[i] = a.poly[i];
		}
		for (; i<=bpwr; i++) {
			add(x.poly[i], a.poly[i], b.poly[i-pwr]);
		}
		for (; i<=a.degree; i++) {
			x.poly[i] = a.poly[i];
		}
	} else {
		d = bpwr;
		//slightly inefficient
		//XXX:find better way
		x = a;
		/*
		for (i=0; i<=b.degree; i++) {
			if (i+pwr > a.degree) {
				x.poly[i+pwr] = b.poly[i];
			} else {
				add(x.poly[i+pwr], x.poly[i+pwr], b.poly[i]);
			}
		}
		*/
		i=pwr;
		while (i <= a.degree) {
			x.poly[i] = b.poly[i-pwr];
			i++;
		}
		while (i <= bpwr) {
			add(x.poly[i], x.poly[i], b.poly[i-pwr]);
			i++;
		}
	}
	for (i=d; i>=0; i--) {
		if (!IsZero(x.poly[i])) break;
	}
	x.degree = i;
}

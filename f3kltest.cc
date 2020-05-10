/* F3kl Test Program
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

#include <iostream.h>

#include <NTL/lzz_pXFactoring.h>

#include "f3kl.h"

void deconvert(zz_pX &f, const F3kl &P) {
	int i;
	clear(f);
	for (i=0; i<=P.degree; i++) {
		if (P.poly[i].b >= 0) {
			SetCoeff(f, i, antilog[P.poly[i].b]);
		}
	}
}

void convert(F3kl &P, const zz_pX &f) {
	int i;
	clear(P);
	for (i=0; i<=deg(f); i++) {
		P.SetCoeff(i, (int) rep(f.rep[i]));
	}
}

int main()
{
	F3kl P, Q, R;
	F3k a;

	f3kl_init(79);
	a.put(7);
	P.SetCoeff(2, a);
	cout << P << endl;
	add(Q, P, P);
	cout << Q << endl;
	mul(R, P, Q);
	cout << R << endl;
	cout << minpoly << endl;

	zz_p::init(3);
	zz_pX f, g, h;
	random(f, 70);
	random(g, 30);
	convert(P, f);
	convert(Q, g);
	zz_pX mp;
	deconvert(mp, minpoly);
	cout << "mp: " << mp << endl;

	F3kl Quo;
	zz_pX hquo;

	cout << "f: " << f << endl;
	cout << "P: " << P << endl;
	cout << "g: " << g << endl;
	cout << "Q: " << Q << endl;
	divrem(Quo, R, P, Q);
	DivRem(hquo, h, f, g);
	cout << "Quo:  " << Quo << endl;
	cout << "hquo: " << hquo << endl;
	cout << "h:    " << h << endl;
	cout << "R:    " << R << endl;
	realmul(R, P, Q);
	mul(h, f, g);
	cout << "h:    " << h << endl;
	cout << "R:    " << R << endl;
	inv(R, P);
	InvMod(h, f, mp);
	cout << "h:    " << h << endl;
	cout << "R:    " << R << endl;
	ZZ n;
	n = 50;
	power(R, P, n);
	PowerMod(h, f, n, mp);
	cout << "h:    " << h << endl;
	cout << "R:    " << R << endl;
	/*
	int i;
	double t0, t1;
	t0 = GetTime();
	for (i=0; i<1000; i++) {
		MulMod(h, f, g, mp);
	}
	t1 = GetTime();
	cout << "NTL: " << (t1 - t0) << endl;

	t0 = GetTime();
	for (i=0; i<1000; i++) {
		mul(R, P, Q);
	}
	t1 = GetTime();
	cout << "special: " << (t1 - t0) << endl;


	zz_pX bigmp;
	BuildIrred(bigmp, 6*79);
	random(f, 6*79);
	random(g, 6*79);
	
	if (0) {
	t0 = GetTime();
	for (i=0; i<1000; i++) {
		MulMod(h, f, g, bigmp);
	}
	t1 = GetTime();
	cout << "NTL big: " << (t1 - t0) << endl;
	}

	ZZ aa, bb ,cc ;
	RandomBits(aa, 11*6*79);
	RandomBits(bb, 11*6*79);
	t0 = GetTime();
	for (i=0; i<1000; i++) {
		mul(cc, aa, bb);
	}
	t1 = GetTime();
	cout << "NTL ZZbig: " << (t1 - t0) << endl;
	*/
	return 0;
}

/*
 * Finds sparse irreducible polynomial of degree l
 */
/*
Copyright (C) 2001 Benjamin Lynn (blynn@cs.stanford.edu)

This file is part of the Stanford short signature system

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

#include <NTL/lzz_pXFactoring.h>

int main()
{
	int larray[] = { 79, 97, 121, 149, 163, 167 };

	zz_p::init(3);
	zz_pX f;

	int li;
	for (li=0; li<6; li++) {

	int i;

	int l = larray[li];

	//BuildIrred(f, l);
	clear(f);
	SetCoeff(f, l, 1);
	SetCoeff(f, 0, 1);

	int count = 1;
	for (;;) {
		int pwr3 = 1;
		i = 1;
		for (;;) {
			int digit = (count % (pwr3 * 3)) / pwr3;
			//cout << digit;
			SetCoeff(f, i, digit);
			pwr3 *= 3;
			if (pwr3 > count) break;
			i++;
		}
		//cout << endl;
		//cout << "Testing " << f << endl;
		if (DetIrredTest(f)) break;
		count++;
	}

	cout << "case " << l << ":\n";
	for (i=0; i<=l; i++) {
		long a = rep(f.rep[i]);
		if (a > 0) {
			cout << "\tminpoly.SetCoeff(" << i << ", " << a << ");" << endl;
		}
	}
	cout << "\tbreak;\n";

	}
	return 0;
}

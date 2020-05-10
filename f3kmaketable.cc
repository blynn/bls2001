/* Precomputes arithmetic lookup tables.
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

#include <NTL/lzz_pXFactoring.h>

int k = 6;
int main()
{
	zz_p::init(3);
	zz_pX f;

	BuildIrred(f, k);

	cout << "//Irred poly = " << f << endl;

	zz_pX g, h;
	int order = 1;
	int i;

	for (i=0; i<k; i++) order *= 3;
	order--;

	bool again;

	//Inefficient, but who cares?
	do {
		again = false;
		random(g, 6);
		for (i=2; i<order; i++) {
			if ((order % i ) == 0) {
				//cout << "factor: " << i << endl;
				PowerMod(h, g, i, f);
				//cout << h << endl;
				if (IsOne(h)) {
					//cout << "not generator\n";
					again = true;
					break;
				}
			}
		}
	} while(again);

	int al[order];
	int dl[order+1];
	zz_pX g2;
	set(g2);
	for (i=0; i<order; i++) {
		int pwr3 = 1;
		int j;
		int c;
		c = 0;
		for (j=0; j<k; j++) {
			c += rep(coeff(g2, j)) * pwr3;
			pwr3 *= 3;
		}
		//cout << c << " = " << g2 << endl;
		al[i] = c;
		dl[c] = i;
		g2 = (g2 * g) % f;
	}
	dl[0] = -order-1;

	int addtable[order];
	for (i=0; i<order; i++) {
		int i1;
		i1 = al[i] + 1;
		if (!(i1 % 3)) i1 -= 3;
		addtable[i] = dl[i1];
	}

	/*
	cout << "check " << g2 <<endl;

	for (i=0; i<order; i++) {
		cout << "DL(" << i << ") = " << dl[i] << endl;
	}
	*/

	int count;
	cout << "int discretelog[] = {";
	count = 0;
	for (i=0; i<order+1; i++) {
		if (!count) cout << endl << "    ";
		count = (count + 1) % 10;
		cout << dl[i];
		if (i < order) {
			cout << ", ";
		}
	}
	cout << endl << "};" << endl;

	cout << endl << "int antilog[] = {";
	count = 0;
	for (i=0; i<order; i++) {
		if (!count) cout << endl << "    ";
		count = (count + 1) % 10;
		cout << al[i];
		if (i < order-1) {
			cout << ", ";
		}
	}
	cout << endl << "};" << endl;

	cout << endl << "int addtable[] = {";
	count = 0;
	for (i=0; i<order; i++) {
		if (!count) cout << endl << "    ";
		count = (count + 1) % 10;
		cout << addtable[i];
		if (i < order-1) {
			cout << ", ";
		}
	}
	cout << endl << "};" << endl;

	return 0;
}

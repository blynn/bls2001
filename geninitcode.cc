/*
Copyright (C) 2001 Benjamin Lynn (blynn@cs.stanford.edu)

This file is part of the short signature system.

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

#include <NTL/ZZ.h>

#include "weil.h"

int main(int argc, char **argv)
{
	int larray[] = { 79, 97, 121, 149, 163, 163, 167 };
	int lsign[] = { -1, 1, 1, 1, 1, -1, 1 };

	int li;
	for (li=0; li<7; li++) {

		cout << "case " << li << ":\n";

		int i;
		int l = larray[li];
		cout << "\tl = " << l << ";" << endl;

		ZZ order, z2;
		power(order, 3, l);
		order++;
		power(z2, 3, (l+1)/2);
		i = l % 12;
		if (i == 1 || i == 11) {
			if (lsign[li] > 0) {
				add(order, order, z2);
			} else {
				sub(order, order, z2);
			}
		} else {
			if (lsign[li] > 0) {
				sub(order, order, z2);
			} else {
				add(order, order, z2);
			}
		}
		
		ZZ pointcount = order;
		//cout << "order = " << order << ";" << endl;
		/*
		for (i=3; i<10000; i+=2) {
			if (!(order % i)) {
				cerr << "found factor " << i << endl;
				order = order / i;
			}
		}
		*/
		PrimeSeq ps;
		long pr = ps.next();
		for (long ii=0; ii<500000; ii++) {
			while (divide(order ,pr)) {
				//cerr << "found factor " << pr << endl;
				order /= pr;
			}
			pr = ps.next();
			/*
			if (pr == 0) {
				cerr << "exhausted!\n";
				exit(1);
			}
			*/
		}
		cout << "\tcurvesign = " << lsign[li] << ";" << endl;
		cout << "\torder = to_ZZ(\"";
		cout << order << "\");" << endl;
		cout << "\tmultiplier = to_ZZ(\"" << pointcount / order << "\");" << endl;
		cout << "\tbreak;\n";

	}
	return 0;
}

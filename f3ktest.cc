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

#include <iostream.h>

#include "f3k.h"

int main()
{
	F3k a, b, c;
	F3k u, rplus, rminus;
	int i;
	for (i=0; i<728; i++) {
		a.b = i;
		sqr(b, a);
		negate(c, b);
		if (IsOne(c)) {
			u = a;
			cout << "u = " << u << endl;
			cout << "i = " << i << endl;
		} else {
			mul(b, b, a);
			add(c, a, a);
			add(b, b, c);
			set(c);
			add(b, b, c);
			if (IsZero(b)) {
				rminus = a;
				cout << "rminus = " << rminus << endl;
				cout << "i = " << i << endl;
			} else {
				add(b, b, c);
				if (IsZero(b)) {
					rplus = a;
					cout << "rplus = " << rplus << endl;
					cout << "i = " << i << endl;
				}
			}
		}
	}

	a.put(14); //6, 1
	b.put(64); //5, 8
	
	add(c, a, b);
	cout << a << " + " << b << " = " << c << endl;
	a.put(7); //1, 2
	b.put(5); //2, 1
	mul(c, a, b);
	cout << a << " * " << b << " = " << c << endl;
	inv(c, a);
	cout << a << " ^-1 = " << c << endl;
	inv(c, c);
	cout << " ^-1 = " << c << endl;

	return 0;
}

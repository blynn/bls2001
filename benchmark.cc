/*
 * Written to produce results for a paper
 * measures how long it takes to verify a signature
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

#include <stdio.h>
#include <NTL/ZZ.h>

#include "weil.h"

F3kl Phi_r;
F3kl Phi_u;

enum {
    use_weil_flag = 0
};

static int randomize_NTL(void)
{
    int n = 20;
    FILE *fp;
    fp = fopen("/dev/urandom", "r");
    if (!fp) {
	return 0;
    }
    char *s;
    s = new char[n];
    if (n != (int) fread(s, 1, n, fp)) {
	delete[] s;
	return 0;
    }
    SetSeed(ZZFromBytes((unsigned char *) s, n));
    delete[] s;
    fclose(fp);
    return 1;
}

void apply_Phi(Point &P)
{
    //P.x = -P.x + Phi_r;
    //P.y = -P.y * Phi_u;

    negate(P.x, P.x);
    add(P.x, P.x, Phi_r);

    negate(P.y, P.y);
    mul(P.y, P.y, Phi_u);
}

int main(int argc, char **argv)
{
    randomize_NTL();
    int choice = 0;

    if (argc > 1) {
	choice = atoi(argv[1]);
    }

    int curvesign;
    int l;
    ZZ order, multiplier;
    switch(choice) {
	case 0:
	    l = 79;
	    curvesign = -1;
	    order = to_ZZ("49269609804781974450852068861184694669");
	    multiplier = to_ZZ("1");
	    break;
	case 1:
	    l = 97;
	    curvesign = 1;
	    order = to_ZZ("2726865189058261010774960798134976187171462721");
	    multiplier = to_ZZ("7");
	    break;
	case 2:
	    l = 121;
	    curvesign = 1;
	    order = to_ZZ("31577918281911659253350468036837879363403758167");
	    multiplier = to_ZZ("170721541921");
	    break;
	case 3:
	    l = 149;
	    curvesign = 1;
	    order = to_ZZ("1159188057595039062484376423422896375894631969050751604761379750689");
	    multiplier = to_ZZ("106393");
	    break;
	case 4:
	    l = 163;
	    curvesign = 1;
	    order = to_ZZ("84268735918094105836318246511533764121140010481130741067443071103148817701717");
	    multiplier = to_ZZ("7");
	    break;
	case 5:
	    l = 163;
	    curvesign = -1;
	    order = to_ZZ("589881151426658740854227725580736348850640632297373414091790995505756623268837");
	    multiplier = to_ZZ("1");
	    break;
	case 6:
	    l = 167;
	    curvesign = 1;
	    order = to_ZZ("6825767609365622572741777967434234893829444442089535769018035530821792872561467");
	    multiplier = to_ZZ("7");
	    break;
	default:
	    cerr << "Bad choice. Must be 0-6\n";
	    exit(1);
	    break;
    }

    cout << "l: " << l << endl;
    cout << "order: " << order << endl;

    weil_init_order(l, order, multiplier);

    F3kl f3kl0, f3kl1, f3kl2;
    clear(f3kl0);
    set(f3kl1);
    add(f3kl2, f3kl1, f3kl1);
    
    F3k temp;

    temp.b = phiuexponent;
    Phi_u.SetCoeff(0, temp);

    if (curvesign > 0) {
	weil_init_curve(f3kl0, f3kl2, f3kl1);
	temp.b = phirplusexponent;
	Phi_r.SetCoeff(0, temp);
    } else {
	weil_init_curve(f3kl0, f3kl2, f3kl2);
	temp.b = phirminusexponent;
	Phi_r.SetCoeff(0, temp);
    }

    int trial_count = 10;
    double total_time = 0;
    for(int i = 0; i < trial_count; i++) {

	Point P, xP;
	ZZ x;
	x = RandomBnd(order);
	random_point(P);

	cout << "random P: " << P << endl;

	cout << "secret x: " << x << endl;
	group_times(xP, x, P);

	cout << "xP: " << xP << endl;

	/*
	F3kl aa, bb;
	sqr(aa, P.x);
	mul(aa, aa, P.x);
	sub(aa, aa, P.x);
	add(aa, aa, f3kl2);
	sqr(bb, P.y);
	sub(aa, aa, bb);
	cout << "check0: " << aa << endl;

	sqr(aa, xP.x);
	mul(aa, aa, xP.x);
	sub(aa, aa, xP.x);
	add(aa, aa, f3kl2);
	sqr(bb, xP.y);
	sub(aa, aa, bb);
	cout << "check0.5: " << aa << endl;

	Point oP;
	group_times(oP, order, P);
	cout << "oP: " << oP << endl;
	*/

	//test: h = H(m)
	Point h;
	random_point(h);
	cout << "H(m) = " << h << endl;

	Point sig;
	group_times(sig, x, h);
	cout << "sig = " << sig << endl;

	F3kl ePS, eHxP;
	/*
	cout <<"check1: ";
	sqr(ePS, sig.x);
	mul(ePS, ePS, sig.x);
	mul(eHxP, f3kl2, sig.x);
	add(ePS, ePS, eHxP);
	add(ePS, ePS, f3kl2);
	sqr(eHxP, sig.y);
	sub(ePS, ePS, eHxP);
	//cout << sig.x * sig.x * sig.x + 2 * sig.x - 1 - sig.y * sig.y << endl;
	cout << ePS << endl;
	*/
	apply_Phi(sig);
	double t0, t1;
	t0 = GetTime();

	if (use_weil_flag) {
	    ePS = weil_pairing(P, sig);
	    //cout << "e(P, Phi(S)) = " << ePS << endl;
	    apply_Phi(xP);
	    eHxP = weil_pairing(h, xP);
	    if (eHxP != ePS) {
		cout << "bug in program!" << endl;
	    }
	} else {
	    ePS = tate_pairing(P, sig);
	    negate(xP.y, xP.y);
	    apply_Phi(xP);
	    mul(eHxP, tate_pairing(h, xP), ePS);
	    if (!IsOne(tate_power(eHxP))) {
		cout << "bug in program!" << endl;
		cout << eHxP << "\n";
	    }
	}
	t1 = GetTime();

	total_time += t1 - t0;
	//cout << "e(H(m), Phi(xP)) = " << eHxP << endl;
	cout << "verification time: " << (t1 - t0) << endl;
    }

    cout << "average verification time: " << total_time / trial_count << endl;

    return 0;
}

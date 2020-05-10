/* Computes Weil pairing, also Tate pairing using Miller's algorithm
 * Modified for short signatures project
 * (Taken from Stanford IBE project)
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

//Notes for myself
//Places that assume x1 = x3 = 0
// - slope_tangent
// - anything to do with vertical lines, e.g. get_vertical
// - group_times?
//I assume a2 = a4 = 0 in the specialized version of group_times
//random_point() no longer used. Instead, we assume p=2 mod 3
//so we can take cube roots easily.

#include <NTL/ZZ.h>
#include "weil.h"

enum {
    //by trial and error, a window size of 4 seems fastest
    windowsize = 4,
    windowsizepower = 7	    //this is 2^windowsize - 1
};

//Fp2 scratchspace
static Fp2* fp2reg0;
static Fp2* fp2reg1;
static Fp2* fp2a;
static Fp2* fp2b;
static Fp2* fp2c;

static int intpower2[windowsize];

static ZZ q; //size of subgroup
static ZZ nonq; //no. of points / size of subgroup
static int maxbitq;
static ZZ p1onq;
static ZZ kpwr;

//coefficients of the elliptic curve
//y^2 = x^3 + a2 * x^2 + a4 * x + a6
static Fp2 *a2, *a4, *a6;

static Fp2 *twicea2;
//precalculate 2 * *a2 for speed

static Fp2 *globalnom, *globaldenom;
//used in miller()

static void slope_tangent(Fp2& res, const Fp2& x, const Fp2& y);
//returns slope of the tangent at (x,y)

static bool miller(Fp2& V, const Point& P, const Point& Phat, const Point& Qhat,
const Point& R1, const Point& R2);

static Point RP1, RP2;

void init_random_point(void)
{
    //pick random R1, R2
    //TODO: precompute these
    //cout << "picking random point 1\n";
    random_point(RP1);
    //cout << "picking random point 2\n";
    random_point(RP2);
}

void weil_init_order(int l, const ZZ& order, const ZZ& multiplier)
{
    f3kl_init(l);
    q = order;
    nonq = multiplier;
    maxbitq = NumBits(q) - 1;

    intpower2[0] = 1;
    int i;
    for (i=1; i<windowsize; i++) intpower2[i] = intpower2[i-1] * 2;

    power(kpwr, 3, 6 * l);
    kpwr--;
    div(kpwr, kpwr, q);
}

void weil_init_curve(Fp2 aa2, Fp2 aa4, Fp2 aa6)
{
    a2 = new Fp2;
    a4 = new Fp2;
    a6 = new Fp2;
    *a2 = aa2;
    twicea2 = new Fp2;
    add(*twicea2, aa2, aa2);
    *a4 = aa4;
    *a6 = aa6;
    if (0) {
	cout << "using curve y^2 = x^3";
	if (!IsZero(aa2)) {
	    cout << " + " << aa2 << "x^2";
	}
	if (!IsZero(aa4)) {
	    cout << " + " << aa4 << "x";
	}
	if (!IsZero(aa6)) {
	    cout << " + " << aa6;
	}
	cout << endl;
    }

    globalnom = new Fp2;
    globaldenom = new Fp2;
    fp2reg0 = new Fp2;
    fp2reg1 = new Fp2;
    fp2a = new Fp2;
    fp2b = new Fp2;
    fp2c = new Fp2;

    init_random_point();
}

Fp2 power(const Fp2& x, const ZZ& n) {
    Fp2 res;
    set(res);

    //use sliding-window method
    Fp2 g[2*windowsizepower+2];
    g[1] = x;
    mul(g[2], x, x);
    int i;

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

    /*
    int m = NumBits(n);
    for(;;) {
	if (bit(n, m)) {
	    mul(res, res, x);
	}
	if (m == 0) break;
	mul(res, res, res);
	m--;
    }
    */
    return res;
}

Fp2 tate_pairing(const Point& P, const Point& Q)
{
    //let Phat = P + R1
    Point Phat, Qhat;
    Fp2 res;

    for(;;) {
	group_plus(Phat, P, RP1);
	group_plus(Qhat, Q, RP2);
	if (!miller(res, P, Phat, Qhat, RP1, RP2)) goto millerfail;
	return res;
millerfail:
	init_random_point();
    }
}

Fp2 tate_power(Fp2& input)
{
    return power(input, kpwr);
}

Fp2 weil_pairing(const Point& P, const Point& Q)
{
    //let Phat = P + R1
    Point Phat, Qhat;
    Fp2 nom, denom;
    Fp2 q;

    for(;;) {
	group_plus(Phat, P, RP1);
	group_plus(Qhat, Q, RP2);
	if (!miller(nom, P, Phat, Qhat, RP1, RP2)) goto millerfail;
	if (!miller(denom, Q, Qhat, Phat, RP2, RP1)) goto millerfail;
	div(q, nom, denom);
	return q;
millerfail:
	init_random_point();

    }
}

bool miller(Fp2& res, const Point& P, const Point& Phat, const Point& Qhat,
	const Point& R1, const Point& R2)
{
    //use sliding-window method
    Fp2 g[2*windowsizepower+2];
    Point bP[2*windowsizepower+2];

    ZZ &n = q;
    Fp2& V = *globalnom;
    Fp2& Vdenom = *globaldenom;
    set(Vdenom);
    set(V);

    //work out f_1
    get_vertical(Qhat, R2, Phat);
    get_line(R2, Qhat, P, R1);
    if (IsZero(V) || IsZero(Vdenom)) return false;
    div(g[1], V, Vdenom);
    bP[1] = P;
    set(Vdenom);
    set(V);

    //work out f_2
    group_plus(bP[2], P, P);
    get_line(Qhat, R2, P, P);
    get_vertical(R2, Qhat, bP[2]);
    if (IsZero(V) || IsZero(Vdenom)) return false;
    //g[2] = g[1] * g[1] * V / Vdenom;
    div(g[2], V, Vdenom);
    mul(g[2], g[2], g[1]);
    mul(g[2], g[2], g[1]);
    set(Vdenom);
    set(V);

    //work out rest of g[], bP[]
    int i;
    for (i=1; i<=windowsizepower; i++) {
	int j = 2 * i + 1;
	int k = j - 2;
	group_plus(bP[j], bP[k], bP[2]);
	get_line(Qhat, R2, bP[k], bP[2]);
	get_vertical(R2, Qhat, bP[j]);
	if (IsZero(V) || IsZero(Vdenom)) return false;
	//g[j] = g[k] * g[2] * V / Vdenom;
	div(g[j], V, Vdenom);
	mul(g[j], g[j], g[k]);
	mul(g[j], g[j], g[2]);
	set(Vdenom);
	set(V);
    }

    Point Z;
    Z.put_O();
    int m;
    m = maxbitq;
    while(m>=0) {
	if (!bit(n, m)) {
	    sqr(V, V);
	    sqr(Vdenom, Vdenom);
	    get_line(Qhat, R2, Z, Z);
	    group_plus(Z, Z, Z);
	    get_vertical(R2, Qhat, Z);
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
	    sqr(V, V);
	    sqr(Vdenom, Vdenom);
	    get_line(Qhat, R2, Z, Z);
	    group_plus(Z, Z, Z);
	    get_vertical(R2, Qhat, Z);
	    for (k=l+1; k<=m; k++) {
		if (bit(n, k)) j += intpower2[k-l];
		sqr(V, V);
		sqr(Vdenom, Vdenom);
		get_line(Qhat, R2, Z, Z);
		group_plus(Z, Z, Z);
		get_vertical(R2, Qhat, Z);
	    }
	    mul(V, V, g[j]);
	    get_line(Qhat, R2, Z, bP[j]);
	    group_plus(Z, Z, bP[j]);
	    get_vertical(R2, Qhat, Z);
	    m = l-1;
	}
    }
    if (IsZero(V) || IsZero(Vdenom)) return false;
    div(res, V, Vdenom);

    /*
    Fp2 check;
    Fp2& f1 = g[1];
    set(Vdenom);
    set(V);
    m = maxbitq;
    for(;;) {
	if (bit(n, m)) {
	    V *= f1;
	    //g
	    get_line(Qhat, R2, Z, P);
	    //h
	    group_plus(Z, Z, P);
	    get_vertical(R2, Qhat, Z);
	}
	if (m == 0) break;
	m--;
	sqr(V, V);
	sqr(Vdenom, Vdenom);
	//g
	get_line(Qhat, R2, Z, Z);
	//h
	group_plus(Z, Z, Z);
	get_vertical(R2, Qhat, Z);
    }
    if (IsZero(V) || IsZero(Vdenom)) return false;
    check = V / Vdenom;

    if (check != res) {
	cerr << "BUG!\n";
	exit(1);
    }
    */
    return true;
}

void random_special(F3kl &x)
{
    int i;
    F3k one, two;
    set(one);
    negate(two, one);
    for (i=0; i<minpolydegree; i++) {
	ZZ c;
	RandomBnd(c, to_ZZ(3));
	if (IsZero(c)) {
	} else if (IsOne(c)) {
	    x.SetCoeff(i, one);
	} else {
	    x.SetCoeff(i, two);
	}
    }
}

ZZ sfs1on4;

void force_subgroup(Point& P)
//make P a point in a subgroup
//achieved by multiplying P appropriately
{
    group_times(P, nonq, P);
}

void random_point(Point& P)
{
    F3kl x, y;

    power(sfs1on4, 3, minpolydegree);
    sfs1on4 = (sfs1on4 + 1) / 4;

    for (;;) {
	random_special(x);
	if (y_from_x(y, x)) break;

    }
    P.set_finite();
    P.put(x, y);
    force_subgroup(P);
}

bool y_from_x(Fp2& y, const Fp2& x)
{
    F3kl x1;

    //x1 = x * (x * (x + *a2) + *a4) + *a6; //Horner's rule
    add(x1, x, *a2);
    mul(x1, x1, x);
    add(x1, x1, *a4);
    mul(x1, x1, x);
    add(x1, x1, *a6);

    //take square root
    power(y, x1, sfs1on4);
    F3kl temp;
    sqr(temp, y);
    sub(temp, temp, x1);
    //must check that it is in fact a square root
    if (IsZero(temp)) {
	return true;
    }
    return false;
}

void Point::random(void)
{
    random_point(*this);
}

void group_times(Point &R, const ZZ &a, const Point &P)
{
    ZZ n;
    bool neg;
    if (a < 0) {
	neg = true;
	n = -a;
    } else {
	neg = false;
	n = a;
    }
    if (P.infinity) {
	R.put_O();
	return;
    }
    if (n == 0) {
	return;
    }
    //use sliding-window method
    Point g[2*windowsizepower+2];
    g[1] = P;
    group_plus(g[2], P, P);
    int i;

    for (i=1; i<=windowsizepower; i++) {
	group_plus(g[2*i+1], g[2*i-1], g[2]);
    }

    int m = NumBits(n) - 1;

    R.put_O();
    while(m>=0) {
	if (!bit(n, m)) {
	    group_plus(R, R, R);
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
	    group_plus(R, R, R);
	    for (k=l+1; k<=m; k++) {
		if (bit(n, k)) j += intpower2[k-l];
		group_plus(R, R, R);
	    }
	    group_plus(R, R, g[j]);
	    m = l-1;
	}
    }
    /*
    check
    Point Rcheck;
    Rcheck.put_O();
    m = NumBits(n) - 1;
    for(;;) {
	if (bit(n, m)) {
	    group_plus(Rcheck, Rcheck, P);
	}
	if (m == 0) break;
	group_plus(Rcheck, Rcheck, Rcheck);
	m--;
    }
    if (R != Rcheck) {
	cerr << "BUG!\n";
	exit(1);
    }
    */

    if (neg) negate(R.y, R.y);
}

void group_plus(Point &R, const Point &P, const Point &Q)
{
    if (P.infinity) {
	R = Q;
	return;
    }
    if (Q.infinity) {
	R = P;
	return;
    }

    R.set_finite();
    //Fp2 lambda, mu;
    Fp2& lambda = *fp2reg0;
    Fp2& mu = *fp2reg1;

    if (P.x == Q.x) {
	negate(lambda, Q.y);
	if (P.y == lambda) {
	    R.put_O();
	} else { //P.y == Q.y
	    //line: Y - (lambda X + mu)

	    slope_tangent(lambda, P.x, P.y);
	    //mu = P.y - lambda * P.x;
	    mul(*fp2a, lambda, P.x);
	    sub(mu, P.y, *fp2a);
	    //R.x = lambda * lambda - P.x - P.x;
	    mul(*fp2a, lambda, lambda);
	    add(*fp2b, P.x, P.x);
	    sub(R.x, *fp2a, *fp2b);
	    //R.y = -lambda * R.x - mu;
	    negate(*fp2a, lambda);
	    mul(*fp2b, *fp2a, R.x);
	    sub(R.y, *fp2b, mu);
	}
    } else {
	//line: Y - (lambda X + mu)

	//lambda = (Q.y - P.y) / (Q.x - P.x);
	sub(*fp2a, Q.y, P.y);
	sub(*fp2b, Q.x, P.x);
	div(lambda, *fp2a, *fp2b);
	//mu = P.y - lambda * P.x;
	mul(*fp2a, lambda, P.x);
	sub(mu, P.y, *fp2a);
	//R.x = lambda * lambda - P.x - Q.x;
	mul(*fp2a, lambda, lambda);
	add(*fp2b, P.x, Q.x);
	sub(R.x, *fp2a, *fp2b);
	//R.y = -lambda * R.x - mu;
	negate(*fp2a, lambda);
	mul(*fp2b, *fp2a, R.x);
	sub(R.y, *fp2b, mu);
    }
}

void get_line(const Point& nom, const Point& denom, const Point& P, const Point& Q)
{
    //cases involving O
    if (P.infinity) {
	get_vertical(nom, denom, Q);
	return;
    }
    if (Q.infinity) {
	get_vertical(nom, denom, P);
	return;
    }
    //check if we need a tangent or vertical
    const Fp2& Px = P.x;
    const Fp2& Py = P.y;
    //Fp2 a, b, c;
    Fp2 &a = *fp2a;
    Fp2 &b = *fp2b;
    Fp2 &c = *fp2c;

    if (Px == Q.x) {
	negate(a, Q.y);
	if (Py == a) {
	    //a = 1; b = 0; c = -P.x;
	    //*globalnom *= nom.x - P.x;
	    sub(a, nom.x, P.x);
	    mul(*globalnom, *globalnom, a);
	    //*globaldenom *= denom.x - P.x;
	    sub(a, denom.x, P.x);
	    mul(*globaldenom, *globaldenom, a);
	    return;
	}
	//it should be
	//a = -slope_tangent(P.x, P.y);
	//b = 1;
	//c = -(P.y + a * P.x);
	//but we multiply by 2*P.y to avoid division

	//a = -Px * (Px + Px + Px + *twicea2) - *a4;
	add(*fp2reg0, Px, Px);
	add(*fp2reg1, Px, *twicea2);
	add(*fp2reg0, *fp2reg0, *fp2reg1);
	mul(*fp2reg0, *fp2reg0, Px);
	add(*fp2reg0, *fp2reg0, *a4);
	negate(a, *fp2reg0);

	//b = Py + Py;
	add(b, Py, Py);

	//c = - b * Py - a * Px;
	mul(*fp2reg0, b, Py);
	mul(*fp2reg1, a, Px);
	add(*fp2reg0, *fp2reg0, *fp2reg1);
	negate(c, *fp2reg0);

	if (nom.infinity) {
	    //*globalnom *= b;
	    mul(*globalnom, *globalnom, b);
	} else {
	    //*globalnom *= a * nom.x + b * nom.y + c;
	    mul(*fp2reg0, a, nom.x);
	    mul(*fp2reg1, b, nom.y);
	    add(*fp2reg0, *fp2reg0, *fp2reg1);
	    add(*fp2reg0, *fp2reg0, c);
	    mul(*globalnom, *globalnom, *fp2reg0);
	}
	if (denom.infinity) {
	    //*globaldenom *= b;
	    mul(*globaldenom, *globaldenom, b);
	} else {
	    //*globaldenom *= a * denom.x + b * denom.y + c;
	    mul(*fp2reg0, a, denom.x);
	    mul(*fp2reg1, b, denom.y);
	    add(*fp2reg0, *fp2reg0, *fp2reg1);
	    add(*fp2reg0, *fp2reg0, c);
	    mul(*globaldenom, *globaldenom, *fp2reg0);
	}
	return;
    }
    //normal simple case

    //it should be
    //a = -(Q.y - P.y) / (Q.x - P.x);
    //b = 1;
    //c = -(P.y + a * P.x);
    //but we'll multiply by Q.x - P.x to avoid division

    //b = Q.x - Px;
    sub(b, Q.x, Px);
    //a = Py - Q.y;
    sub(a, Py, Q.y);
    //c = - b * Py - a * Px;
    mul(c, b, Py);
    mul(*fp2reg0, a, Px);
    add(c, c, *fp2reg0);
    negate(c, c);

    if (nom.infinity) {
	//*globalnom *= b;
	mul(*globalnom, *globalnom, b);
    } else {
	//*globalnom *= a * nom.x + b * nom.y + c;
	mul(*fp2reg0, a, nom.x);
	mul(*fp2reg1, b, nom.y);
	add(*fp2reg0, *fp2reg0, *fp2reg1);
	add(*fp2reg0, *fp2reg0, c);
	mul(*globalnom, *globalnom, *fp2reg0);
    }
    if (denom.infinity) {
	//*globaldenom *= b;
	mul(*globaldenom, *globaldenom, b);
    } else {
	//*globaldenom *= a * denom.x + b * denom.y + c;
	mul(*fp2reg0, a, denom.x);
	mul(*fp2reg1, b, denom.y);
	add(*fp2reg0, *fp2reg0, *fp2reg1);
	add(*fp2reg0, *fp2reg0, c);
	mul(*globaldenom, *globaldenom, *fp2reg0);
    }
}

void get_vertical(const Point& nom, const Point& denom, const Point& P)
{
    if (nom.infinity || denom.infinity) {
	clear(*globalnom);
	clear(*globaldenom);
	return;
    }
    if (P.infinity) {
	//a = b = 0; c = 1;
	return;
    }
    //a = 1; b = 0; c = -P.x;
    //*globalnom *= nom.x - P.x;
    sub(*fp2reg0, nom.x, P.x);
    mul(*globalnom, *globalnom, *fp2reg0);
    //*globaldenom *= denom.x - P.x;
    sub(*fp2reg0, denom.x, P.x);
    mul(*globaldenom, *globaldenom, *fp2reg0);
}

void slope_tangent(Fp2& res, const Fp2& x, const Fp2& y)
{
    //res = (x * (x + x + x + *twicea2) + *a4) / (y + y);
    /*
    add(*fp2reg0, x, x);
    add(*fp2reg1, *twicea2, x);
    add(*fp2reg0, *fp2reg0, *fp2reg1);
    mul(*fp2reg0, x, *fp2reg0);
    add(*fp2reg0, *fp2reg0, *a4);
    add(*fp2reg1, y, y);
    div(res, *fp2reg0, *fp2reg1);
    */
    //char 3 version:
    div(res, *a4, y);
    add(res, res, res);
}

long operator !=(const Point& P, const Point& Q)
{
    return (P.infinity != Q.infinity || P.x != Q.x || P.y != Q.y);
}
long operator ==(const Point& P, const Point& Q)
{
    return (P.infinity == Q.infinity && P.x == Q.x && P.y == Q.y);
}

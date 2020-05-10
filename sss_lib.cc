/*
 * Stanford short signature library
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
#include <openssl/evp.h>

#include "weil.h"

static ZZ pwr3l;

unsigned char *workchar;

F3kl Phi_r;
F3kl Phi_u;

static EVP_MD_CTX mdctx;
const static EVP_MD *md;

//base-36
/*
const static int digit[] = {
	'0', '1', '2', '3', '4', '5', '6', '7', '8', '9',
	'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J',
	'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T',
	'U', 'V', 'W', 'X', 'Y', 'Z' };
const static int base = 36;
*/

//base-64
const static int digit[] = { 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J',
	'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T',
	'U', 'V', 'W', 'X', 'Y', 'Z',
	'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j',
	'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't',
	'u', 'v', 'w', 'x', 'y', 'z',
	'0', '1', '2', '3', '4', '5', '6', '7', '8', '9',
	'+', '/' };
const static int base = 64;

static int digitinverse[128];

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

ZZ order;

void SSS_init(int choice)
{
    OpenSSL_add_all_algorithms();
    md = EVP_get_digestbyname("SHA1");
    randomize_NTL();

    int curvesign;
    int l;
    ZZ multiplier;
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

    power(pwr3l, 3, l);
    //cout << "l: " << l << endl;
    //cout << "order: " << order << endl;
    workchar = new unsigned char[NumBytes(pwr3l) + 5];

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

    int i;
    for(i=0; i<128; i++) {
	digitinverse[i] = -1;
    }
    for(i=0; i<base; i++) {
	digitinverse[digit[i]] = i;
    }
    digitinverse[0] = base;
}

void number_from_hash(ZZ& x, const ZZ& limit, const unsigned char *c, const int len)
{
    ZZ z;
    long countlen;
    ZZ count;
    set(count);
    long i = 0;
    long bits = NumBits(limit);
    ZZFromBytes(z, c, len);
    long zbits = NumBits(z);
    //cout << "bits: " << bits << endl;
    x = 0;
    do {
	int j;
	for (j=0; j<zbits; j++) {
	    if (bit(z, j)) {
		SetBit(x, i);
	    }
	    i++;
	}
	bits -= zbits;
	//cout << "bits: " << bits << endl;
	if (bits <= 0) break;

	countlen = NumBits(count);
	for (j=0; j<countlen; j++) {
	    if (bit(count, j)) {
		SetBit(x, i);
	    }
	    i++;
	}
	bits -= countlen;
	//cout << "bits: " << bits << endl;
	count++;
    } while (bits > 0);
    //cout << "hash: " << x << endl;
    while (x > limit) {
	//cout << x << " is too big (" << limit << ")\n";
	SwitchBit(x, NumBits(x) - 1);
	//cout << "new hash: " << x << endl;
    }
}

void hash_to_point(Point &P, const unsigned char *c, const int len)
{
    int n = 20;
    unsigned char *md_value;
    md_value = new unsigned char[n];
    unsigned int md_len;
    unsigned char *countchar;
    int countlen;
    countlen = 10;
    countchar = new unsigned char[countlen];
    ZZ count;
    clear(count);

    F3k one, two;
    set(one);
    negate(two, one);

    F3kl x, y;
    for(;;) {

	EVP_DigestInit(&mdctx, md);
	EVP_DigestUpdate(&mdctx, c, len);
	BytesFromZZ(countchar, count, countlen);
	EVP_DigestUpdate(&mdctx, countchar, countlen);
	EVP_DigestFinal(&mdctx, md_value, &md_len);

	ZZ r;
	number_from_hash(r, pwr3l, md_value, md_len);
	//cout << "r: " << r << "\n";

	clear(x);
	int i;
	int m;
	for (i=0; i<minpolydegree; i++) {
	    m = r % 3;
	    if (m == 1) {
		x.SetCoeff(i, one);
	    } else if (m == 2) {
		x.SetCoeff(i, two);
	    }
	    r /= 3;
	}
	count++;

	if (y_from_x(y, x)) break;
    }

    delete[] md_value;
    delete[] countchar;
    P.put(x, y);
    force_subgroup(P);
}

//XXX: order of arguments differs from rest of program
//I normally have thing(s) I'm assigning to first.
void string_from_ZZ(const ZZ& z, unsigned char*& data, int& len)
{
    int n;
    n = NumBytes(z);
    if ((n >> 16) > 0) {
	cerr << "z takes more than 16 bytes to represent\n";
	exit(1);
    }
    len = n + 2;
    data = new unsigned char[len];
    data[0] = n >> 8;
    data[1] = n & 255;
    BytesFromZZ(&data[2], z, n);
}

int ZZ_from_string(ZZ& x, const unsigned char *xtext, const int xlen)
{
    int i, n;
    if (xlen < 2) {
	cerr << "bad serialization of integer" << endl;
	return 0;
    }
    n = (xtext[0] << 8) + xtext[1];
    i = 2;
    if (n + 2 > xlen) {
	cerr << "bad serialization of integer" << endl;
	return 0;
    }
    int j;
    for (j=0; j<n; j++) {
	workchar[j] = xtext[i];
	i++;
    }
    ZZFromBytes(x, workchar, n);
    return i;
}

void ZZFromF3kl(ZZ& z, const F3kl &f)
{
    ZZ pwr3;
    set(pwr3);
    int i;
    clear(z);
    for (i=0; i<=f.degree; i++) {
	int temp;
	temp = f.poly[i].b;
	if (temp < 0) { //i.e. zero
	} else if (temp == 0) { //i.e. one
	    z += pwr3;
	} else {
	    z += 2 * pwr3;
	}
	pwr3 *= 3;
    }
}

void F3klFromZZ(F3kl &f, const ZZ& z1)
{
    ZZ z = z1;
    int i;
    i = 0;
    F3k one, two;
    set(one);
    negate(two, one);
    while(!IsZero(z) && i < minpolydegree) {
	int m;
	m = z % 3;
	if (m == 1) {
	    f.SetCoeff(i, one);
	} else if (m == 2) {
	    f.SetCoeff(i, two);
	}
	z /= 3;
	i++;
    }
}

void encode_point(char*&s, const Point& P)
{
    int i, l;
    ZZ x;

    RandomBits(x, 320);
    s = new char[2 * NumBytes(x) + 1];
    l = 0;
    while (!IsZero(x)) {
	i = x % base;
	s[l] = digit[i];
	x = x / base;
	l++;
    }
    s[l] = '\0';
    cout << "320: " << s << endl;

    ZZFromF3kl(x, P.x);
    s = new char[2 * NumBytes(x) + 1];
    l = 0;
    while (!IsZero(x)) {
	i = x % base;
	s[l] = digit[i];
	x = x / base;
	l++;
    }
    s[l] = '\0';
}

int decode_point(Point &P, const char *c)
{
    ZZ x, y;
    ZZ pwrbase;
    int i = 0;
    set(pwrbase);
    clear(x);
    for(;;) {
	if (c[i] < 0) {
	    return 0;
	} else {
	    int j = digitinverse[c[i]];
	    if (j == base) {
		i++;
		break;
	    } else if (j < 0) {
		return 0;
	    } else {
		x += pwrbase * j;
		pwrbase *= base;
		i++;
	    }
	}
    }
    set(pwrbase);
    /*
    clear(y);
    for(;;) {
	if (c[i] <= '9' && c[i] >= '0') {
	    y += pwrbase * (c[i] - '0');
	} else if (c[i] >= 'A' && c[i] <= 'Z') {
	    y += pwrbase * (c[i] - 'A' + 10);
	} else if (!c[i]) {
	    break;
	} else {
	    return 0;
	}
	pwrbase *= 36;
	i++;
    }
    F3klFromZZ(P.y, y);
    */
    F3klFromZZ(P.x, x);
    if (!y_from_x(P.y, P.x)) {
	return 0;
    }
    return i;
}

void old_encode_point(unsigned char*& c, int &len, const Point &P)
{
    unsigned char *s1, *s2;
    int n1, n2;
    ZZ z;

    ZZFromF3kl(z, P.x);
    string_from_ZZ(z, s1, n1);
    ZZFromF3kl(z, P.y);
    string_from_ZZ(z, s2, n2);
    len = n1 + n2;
    c = new unsigned char[len];
    memcpy(c, s1, n1);
    memcpy(&c[n1], s2, n2);
    delete[] s1, s2;
}

int old_decode_point(Point &P, const unsigned char *c, int len)
{
    if (len < 2) {
	cerr << "bad serialization of integer" << endl;
	return 0;
    }
    int i, n;
    n = (c[0] << 8) + c[1];
    i = 2;

    if (n + 2 > len) {
	cerr << "bad serialization of integer" << endl;
	return 0;
    }
    int j;
    for (j=0; j<n; j++) {
	workchar[j] = c[i];
	i++;
    }
    ZZ x;
    ZZFromBytes(x, workchar, n);
    n = (c[i] << 8) + c[i + 1];
    i += 2;
    for (j=0; j<n; j++) {
	workchar[j] = c[i];
	i++;
    }
    ZZ y;
    ZZFromBytes(y, workchar, n);
    F3klFromZZ(P.x, x);
    F3klFromZZ(P.y, y);

    return i;
}

void BytesFromPK(unsigned char*& pk, int& pklen, Point P, Point xP)
{
    unsigned char *c1, *c2;
    int l1, l2;
    ZZ z;

    old_encode_point(c1, l1, P);
    old_encode_point(c2, l2, xP);

    pklen = l1 + l2;
    pk = new unsigned char[pklen];
    memcpy(pk, c1, l1);
    memcpy(&pk[l1], c2, l2);
    delete[] c1, c2;
}

int PKFromBytes(Point& P, Point& xP, const unsigned char *pk, int pklen)
{
    int l1, l2;
    l1 = old_decode_point(P, pk, pklen);
    if (l1 == 0 || l1 >= pklen) return 0;
    l2 = old_decode_point(xP, &pk[l1], pklen - l1);
    if (l2 == 0) return 0;
    return l1 + l2;
}

void SSS_generate_keypair(unsigned char*& pk, int& pklen,
	unsigned char*& sk, int& sklen)
{
    ZZ x;
    Point P, xP;
    x = RandomBnd(order);
    random_point(P);

    //cout << "rp P: " << P << endl;
    //cout << "x: " << x << endl;
    group_times(xP, x, P);

    //cout << "xP: " << xP << endl;
    string_from_ZZ(x, sk, sklen);
    BytesFromPK(pk, pklen, P, xP);
}

//void SSS_sign(unsigned char*& sig, int& siglen, const unsigned char *c, int len,
	//const unsigned char *sk, int sklen)
void SSS_sign(char*& sig, const unsigned char *c, int len,
	const unsigned char *sk, int sklen)
{
    ZZ x;
    Point Sig;
    hash_to_point(Sig, c, len);
    ZZ_from_string(x, sk, sklen);
    group_times(Sig, x, Sig);
    encode_point(sig, Sig);
}

//bool SSS_verify(const unsigned char* sig, int siglen, const unsigned char *c,
	//int len, const unsigned char *pk, int pklen);
bool SSS_verify(const char* sig, const unsigned char *c,
	int len, const unsigned char *pk, int pklen)
{
    Point H;
    hash_to_point(H, c, len);
    Point Sig;
    if (!decode_point(Sig, sig)) return false;
    apply_Phi(Sig);
    Point P, xP;
    F3kl eHxP, ePS;
    if (!PKFromBytes(P, xP, pk, pklen)) return false;
    //double t0, t1;
    //t0 = GetTime();
    ePS = tate_power(tate_pairing(P, Sig));
    apply_Phi(xP);
    eHxP = tate_power(tate_pairing(H, xP));
    //t1 = GetTime();
    //cout << "verification time: " << (t1 - t0) << " seconds" << endl;
    if (eHxP != ePS) {
	F3kl prod;
	mul(prod, eHxP, ePS);
	if (!IsOne(prod)) {
	    return false;
	}
    }
    return true;
}

int SSS_base36convert(char*&s, int&l, unsigned char* c, int len)
{
    if (len < 2) {
	cerr << "bad serialization of integer" << endl;
	return 0;
    }
    int i, n;
    n = (c[0] << 8) + c[1];
    i = 2;

    if (n + 2 > len) {
	cerr << "bad serialization of integer" << endl;
	return 0;
    }
    int j;
    for (j=0; j<n; j++) {
	workchar[j] = c[i];
	i++;
    }
    ZZ x;
    ZZFromBytes(x, workchar, n);
    n = (c[i] << 8) + c[i + 1];
    i += 2;
    for (j=0; j<n; j++) {
	workchar[j] = c[i];
	i++;
    }
    ZZ y;
    ZZFromBytes(y, workchar, n);

    s = new char[2 * (NumBytes(y) + NumBytes(x) + 1)];
    l = 0;
    while (!IsZero(x)) {
	i = x % 36;
	if (i<10) {
	    s[l] = i + '0';
	    cout << i;
	} else {
	    s[l] = 'A' + (i-10);
	    cout << (char) ('A' + (i-10));
	}
	x = x / 36;
	l++;
    }
    /*
    cout << "-";
    s[l] = '-';
    l++;
    while (!IsZero(y)) {
	i = y % 36;
	if (i<10) {
	    s[l] = i + '0';
	    cout << i;
	} else {
	    s[l] = 'A' + (i-10);
	    cout << (char) ('A' + (i-10));
	}
	y = y / 36;
	l++;
    }
    */
    cout << endl;
    s[l] = '\0';

    return i;
}

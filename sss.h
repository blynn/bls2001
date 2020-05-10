/*
 * Stanford short signature library header file
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

void SSS_init(int choice);
void SSS_generate_keypair(unsigned char*& pk, int& pklen,
		unsigned char*& sk, int& sklen);
void SSS_sign(char*& sig, const unsigned char *c, int len,
		const unsigned char *sk, int sklen);
//void SSS_sign(unsigned char*& sig, int& siglen, const unsigned char *c, int len,
		//const unsigned char *sk, int sklen);
bool SSS_verify(const char* sig, const unsigned char *c,
		int len, const unsigned char *pk, int pklen);
//bool SSS_verify(const unsigned char* sig, int siglen, const unsigned char *c,
		//int len, const unsigned char *pk, int pklen);

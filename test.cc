#include <iostream.h>
#include <iomanip.h>
#include <stdlib.h> //for atoi()
#include <NTL/tools.h>
#include "sss.h"

int main(int argc, char **argv) {
    int choice;
    choice = 0;
    if (argc > 1) {
	choice = atoi(argv[1]);
	if (choice < 0 || choice > 6) {
	    cerr << "bad choice: must be 0-6\n";
	    exit(1);
	}
    }
    cout << "choice: " << choice << endl;

    SSS_init(choice);
    unsigned char *pk, *sk;
    int pklen, sklen;

    SSS_generate_keypair(pk, pklen, sk, sklen);

    char *sig;
    SSS_sign(sig, (unsigned char *) "Ben Lynn", 8, sk, sklen);
    cout << sig << endl;

    if (SSS_verify(sig, (unsigned char *) "Ben Lynn", 8, pk, pklen)) {
	cout << "sig verified!\n";
    }
    return 0;
}

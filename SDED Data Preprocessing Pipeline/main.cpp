#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <ctime>

#include "SUPREME_PPP.h"

int main(int argc, char **argv)
{
	srand((unsigned int)time(NULL));
	for (int cirnum = 0; cirnum < 1; cirnum++)
	{
		if (argc > 1)
			SUPREME_PPP_Processing(argv[1]);
		else
			SUPREME_PPP_Processing(NULL);
	}
	return 1;

}

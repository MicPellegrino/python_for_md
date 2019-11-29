#!/usr/bin/awk -f

BEGIN {
    # Defaults
    ni = 500;
    nj = 20;
    nk = 10;
    sp = 2.7;

    if (ARGV[1] == "help") {print "usage: make_lj ni nj nk sp"; exit}
    if (ARGC >= 2) ni = ARGV[1]
    if (ARGC >= 3) nj = ARGV[2]
    if (ARGC >= 4) nk = ARGV[3]
    if (ARGC >= 5) sp = ARGV[4]

    dx = sp;
    dy = sp*sqrt(3/4);
    dz = sp*sqrt(2/3);

    dx_y = dx/2;
    dx_z = dx/2;
    dy_z = dy/3;

    x0 = dx/4;
    y0 = dy/6;
    z0 = dz/2;

    printf("CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1\n",
           ni*dx,nj*dy,nk*dz,90.0,90.0,90.0);

    n = 1;
    for(k=0; k<nk; k++) {
	if (k < nk - 1)
	    a = "CUB";
	else
	    a = "CUS";
	z = z0 + k*dz;
	for(i=0; i<ni; i++) {
	    for(j=0; j<nj; j++) {
		y = y0 + j*dy + k*dy_z;
		x = x0 + i*dx + j*dx_y + k*dx_z;
		printf("%-6s%5d  %-3s%4s %5d    %8.3f%8.3f%8.3f\n",
		       "ATOM",n % 100000,a,"SUB",1,x,y,z);
		n++;
	    }
	}
    }
}

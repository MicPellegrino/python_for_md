#!/usr/bin/awk -f

BEGIN {
    # number of molecules in x, y, z respectively
    ni = 220;
    nj = 12;
    nk = 1;

    if (ARGV[1] == "help") {print "usage: make_quad ni nj nk"; exit}
    if (ARGC >= 2) ni = ARGV[1]
    if (ARGC >= 3) nj = ARGV[2]
    if (ARGC >= 4) nk = ARGV[3]

    # names
    as = "SI";
    ao1 = "O1";
    ao2 = "O2";

    # spacing constant
    sp = 4.50;

    d_so = 1.51;

    # calculate spacing
    dx = sp;
    dy = sp*sqrt(3/4);
    dz = sp*sqrt(2/3);

    # row spacing
    dx_y = dx/2;
    #dx_z = dx/2;
    #dy_z = dy/3;
    dx_z = 0;
    dy_z = 0;

    # starting positions
    x0 = dx/4;
    y0 = dy/6;
    #z0 = dz/2;
    z0 = 3;

    # print system information
    printf("CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1\n",
           ni*dx,nj*dy,nk*dz,90.0,90.0,90.0);

    n = 1;
    # loop over all layers
    for (k = 0; k < nk; k++)
    {
        # check if bulk or not
        if (k < nk - 1)
        {
            a = "CUB";
        }
        else
        {
            a = "CUS";
        }

        z = z0 + k*dz;
        # loop over rows and columns
        for(i = 0; i < ni; i++)
        {
            for(j = 0; j < nj; j++)
            {
                y = y0 + j*dy + k*dy_z;
                x = x0 + i*dx + j*dx_y + k*dx_z;
                printf("%-6s%5d  %-3s%4s %5d    %8.3f%8.3f%8.3f\n",
                       "ATOM", n % 100000, ao1, "SUB", 1, x, y, z+d_so);
                printf("%-6s%5d  %-3s%4s %5d    %8.3f%8.3f%8.3f\n",
                       "ATOM", n % 100000, as, "SUB", 1, x, y, z);
                printf("%-6s%5d  %-3s%4s %5d    %8.3f%8.3f%8.3f\n",
                       "ATOM", n % 100000, ao2, "SUB", 1, x, y, z-d_so);
                n++;
            }
        }
    }
}

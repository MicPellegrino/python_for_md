#!/usr/bin/awk -f

BEGIN {
    if (ARGC == 3) {
        cut2 = ARGV[2]*ARGV[2];
        ARGV[2] = "";
    } else {
        cut2 = 0;
    }

    nat_per_mol = 3;

    na = 0;
}

{
    if ($1 == "CRYST1") {
	cx = $2/2;
	cy = 0;
	cz = $4/2;
	if (cut2 == 0) {
	    # cut2 = cx*cx;
      cut2 = (cz*cz)/2;
	}
    }
    if ($1=="ATOM") {
	aline[na] = $0;
	na++;
	if (na == 1)
	{
	    #x = $6;
	    #y = $7;
	    #z = $8;
	    x = substr($0,31,8);
	    y = substr($0,39,8);
	    z = substr($0,47,8);
	    r2 = (x-cx)*(x-cx) + (z-cz)*(z-cz);
	}
	if (na == nat_per_mol) {
	    if (r2 < cut2) {
		for(a=0; a<na; a++) {
		    print aline[a];
		}
	    }
	    na = 0;
	}
    } else {
	print $0;
    }
}

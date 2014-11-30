# Awk program called in utest.sh
{
    # v = version 0 for min, 1 for max
    # t = target lambda
    # l = current lambda
    # r = percent difference between l and t
    # e = epsilon or comparison threshold on r
    l = $3
    r = (l - t) / t;
    if (r < 0) {
        r = -r; # abs()
    }    
    if (v == 0) {
        v = "max"
    } else {
        v = "min"
    }
    s = "program=" p " ver=" v " val=" l " target=" t " %diff=" r " threshold=" e;
    if (r <= e) {
        print "Passed:", s;
    } else {
        print "Failed:", s;
    }
}

# EOF

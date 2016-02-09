while True:
    conversion = raw_input("conversion? (0.0-1.0) ")
    try:
        conversion = float(conversion)
    except:
        print "invalid conversion"        

    if 0.0 < conversion < 1.0:
        break
    else:
        print "out of range"

while True:
    r_major = raw_input("r_major, the amount of the major isotope (>0.0)? ")
    try:
        r_major = float(r_major)
    except:
        print "invalid r_major"     
        continue

    if r_major < 0.0:
        print "out of range"
        continue
    
    r_minor = raw_input("r_minor, the amount of the minor isotope (>0.0)? ")
    try:
        r_minor = float(r_minor)
    except:
        print "invalid r_minor"       
        continue

    if r_minor < 0.0:
        print "out of range"
        continue
 
    if r_minor > r_major:
        print "major must be greater than minor"
        continue

    break

while True:
    KIE = raw_input("KIE? ")
    try:
        KIE = float(KIE)
    except:
        print "invalid conversion"        

    if KIE > 0.0:
        break
    else:
        print "out of range"

r0 = r_major/r_minor
print "starting isotopic ratio is %.5f" % r0

print "\n==recovered sm=="
enhancement = (1.0 - conversion) ** (KIE-1)
print "isotopic ratio changes by a factor of %.5f (%.5f)" % (enhancement, 1.0/enhancement)

r = r0 * enhancement
print "isotopic ratio is now %.5f" % r

major = (r * 100.0) / (r + 1.0)
minor = 100.0-major

print "major isotope: %.5f%%" % major
print "minor isotope: %.5f%%" % minor

print "\n==recovered pdt=="
enhancement = (1.0 - conversion) ** ( 1.0 / KIE )
enhancement = (1.0 - enhancement) / conversion
print "isotopic ratio changes by a factor of %.5f (%.5f)" % (enhancement, 1.0/enhancement)
r = r0 * enhancement
print "isotopic ratio is now %.5f" % r

major = (r * 100.0) / (r + 1.0)
minor = 100.0-major

print "major isotope: %.5f%%" % major
print "minor isotope: %.5f%%" % minor



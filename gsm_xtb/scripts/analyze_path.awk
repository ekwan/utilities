# usage: awk -f analyze_path.awk stringfile.xyz0001

{
    if ( NF == 1 ) {
        getline
        n_geoms++
        n_atoms=0
        energy[n_geoms] = $1
    }
    else {
        n_atoms++
        x[n_geoms,n_atoms]=$2
        y[n_geoms,n_atoms]=$3
        z[n_geoms,n_atoms]=$4
    }
}

END {
    print "node   energy   dist1   dist2"
    print "-----------------------------"
    for (i=1; i <= n_geoms; i++) {
        #distance1 = distance(i,6,1)
        distance1 = dihedral(i,1,5,8,11)
        #distance2 = distance(i,1,5)
        printf " %2d    %5.2f   %6.3f  %6.3f\n", i, energy[i], distance1, distance2
    }
}

function distance(i,atom1,atom2,   temp,x1,x2,y1,y2,z1,z2) {
    x1=x[i,atom1]
    y1=y[i,atom1]
    z1=z[i,atom1]
    x2=x[i,atom2]
    y2=y[i,atom2]
    z2=z[i,atom2]
    temp=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2)
    temp=sqrt(temp)
    return temp
}

function angle(i,atomNumber1,atomNumber2,atomNumber3, value) {
    value=((-distance(i,atomNumber1,atomNumber3)^2+distance(i,atomNumber1,atomNumber2)^2+distance(i,atomNumber2,atomNumber3)^2)/(2*distance(i,atomNumber1,atomNumber2)*distance(i,atomNumber2,atomNumber3)))
   return acos(value)
}

function asin(x) { return (180/3.141592)*atan2(x, sqrt(1-x*x)) }

function acos(x) { return (180/3.141592)*atan2(sqrt(1-x*x), x) }

function atan(x) { return (180/3.141592)*atan2(x,1) }

function dihedral(i,atomNumber1,atomNumber2,atomNumber3,atomNumber4) {
    x1=x[i,atomNumber1]
    y1=y[i,atomNumber1]
    z1=z[i,atomNumber1]
    x2=x[i,atomNumber2]
    y2=y[i,atomNumber2]
    z2=z[i,atomNumber2]
    x3=x[i,atomNumber3]
    y3=y[i,atomNumber3]
    z3=z[i,atomNumber3]
    x4=x[i,atomNumber4]
    y4=y[i,atomNumber4]
    z4=z[i,atomNumber4]

    B1x=x2-x1
    B1y=y2-y1
    B1z=z2-z1

    B2x=x3-x2
    B2y=y3-y2
    B2z=z3-z2

    B3x=x4-x3
    B3y=y4-y3
    B3z=z4-z3

   modB2=sqrt((B2x^2)+(B2y^2)+(B2z^2))
# yAx is x-coord. etc of modulus of B2 times B1
   yAx=modB2*(B1x)
   yAy=modB2*(B1y)
   yAz=modB2*(B1z)
# CP2 is the crossproduct of B2 and B3
   CP2x=(B2y*B3z)-(B2z*B3y)
   CP2y=(B2z*B3x)-(B2x*B3z)
   CP2z=(B2x*B3y)-(B2y*B3x)
   termY=((yAx*CP2x)+(yAy*CP2y)+(yAz*CP2z))
# CP is the crossproduct of B1 and B2
   CPx=(B1y*B2z)-(B1z*B2y)
   CPy=(B1z*B2x)-(B1x*B2z)
   CPz=(B1x*B2y)-(B1y*B2x)
   termX=((CPx*CP2x)+(CPy*CP2y)+(CPz*CP2z))
   dihed4=(180/3.141592)*atan2(termY,termX)
   return dihed4
}


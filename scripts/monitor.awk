# monitors a single g09 opt job
# don't use for scans or anything else
# awk -f monitor.awk job.out

BEGIN { iterations = 1
scf_start = -1
scf_end = -1
gradient_end = -1
}

NR != FNR {
    print "only supports one output file at a time, ignoring additional input files"
    exit
}

/Standard orientation/,/Rotational constants/ {
    if (match($0,"[a-zA-Z]") == 0 && NF==6)
        {
            atomSymbol[$1]=$2
            x[iterations,$1]=$4
            y[iterations,$1]=$5
            z[iterations,$1]=$6
            numberOfAtoms=$1
        }
}

/SCF Done/ {
    energy[iterations]=$5
    if ( $5 < minEnergy || minEnergy == 0.0 )
        minEnergy = $5
    iterations++
}

/ Maximum Force/ && NF == 5 {
    maxForces[iterations-1]=$3
    getline
    RMSforces[iterations-1]=$3
    getline
    maxDisplacements[iterations-1]=$3
    getline
    RMSdisplacements[iterations-1]=$3
}

/Normal termination/ {
    finished++
}

/^ Frequencies --/ {
    didFrequencies=1
    for (i=3; i <= NF; i++)
        {
            if ($i < 0)
                {
                    numberOfImaginaries++
                    imaginaries[numberOfImaginaries]=$i
                }
        }
}

# start of an scf cycle
/Leave Link  202/ {
    scf_start = substr($0, 29, 11)
}

# end of an scf cycle
/Leave Link  502/ {
    scf_end = substr($0, 29, 11)
}

# end of a gradient cycle
/Leave Link  103/ {
    gradient_end = substr($0, 29, 11)
    if (iterations == 1)
        next
    scf_time[iterations-1] = time_delta(scf_start, scf_end)
    gradient_time[iterations-1] = time_delta(scf_end, gradient_end)
    scf_start = -1
    scf_end = -1
    gradient_end = -1
}

END {
    print "iterations  rel_energy   rms_force    max_force   rms_disp   max_disp   scf_time   grad_time   bond1    bond2"
    for (i=1; i < iterations; i+=1)
        {
            relEnergy = (energy[i] - minEnergy)*627.509469
            maxForce = i in maxForces ? sprintf("%9.1f", maxForces[i]*1E5) : ""
            RMSforce = i in RMSforces ? sprintf("%9.1f", RMSforces[i]*1E5) : ""
            maxDisplacement = i in maxDisplacements ? sprintf("%9.1f", maxDisplacements[i]*1E5): ""
            RMSdisplacement = i in RMSdisplacements ? sprintf("%9.1f", RMSdisplacements[i]*1E5): ""
            bond1 = distance(i, 6,9)
            bond2 = distance(i, 9,12)
            printf "  %3d       %9.4f  %9s    %9s   %9s   %9s    %5.2f      %5.2f     %6.3f   %6.3f\n", i, relEnergy, RMSforce, maxForce, RMSdisplacement, maxDisplacement, scf_time[i], gradient_time[i], bond1, bond2
        }
    printf "Normal terminations: %d\n", finished
    if ( didFrequencies == 1 )
        {
            for (i=1; i <= numberOfImaginaries; i++)
                imaginaryString += imaginaries[i] " "
            if (numberOfImaginaries>0)
                printf "Imaginaries: %s\n", imaginaryString
            else
                print "All positive frequencies."
        }
    else
        print "No frequencies calculated."
}

# difference between two times in minutes
# time1 = start, time2 = end
function time_delta(time1, time2,       fields1, fields2) {
    split(time1, fields1, /[ :]/)
    split(time2, fields2, /[ :]/)
    # if we happen to go across months
    if (fields1[1] > fields2[1])
        fields2[1] = fields1[1]+1
    days = (fields2[1] - fields1[1])*24*60
    hours = (fields2[2] - fields1[2])*60
    minutes = fields2[3] - fields1[3]
    seconds = (fields2[4] - fields1[4])/60
    #print time1, time2, days+hours+minutes+seconds
    return days+hours+minutes+seconds
}

function calculateRMS(i,j,        k,RMS,deltaX,deltaY,deltaZ) {
    if (i<1 || j<1)
        return 0.0
    for (k=1; k <= numberOfAtoms; k++)
        {
            deltaX = (x[i,k]-x[j,k])^2
            deltaY = (y[i,k]-y[j,k])^2
            deltaZ = (y[i,k]-z[j,k])^2
            RMS += deltaX + deltaY + deltaZ
        }
    RMS = sqrt(RMS)/numberOfAtoms
    return RMS
}

# finds out if a torsion angle is within tolerance of a target
# angles are assumed be between [-180,180]
function near(value,target,tolerance,      temp1, temp2, diff1, diff2)
{
    temp1 = value+180.0
    temp2 = target+180.0
    diff1 = abs(temp1-temp2)
    diff2 = 360-diff1
    if ( min(diff1,diff2) <= tolerance )
        return 1 # true
    return 0     # false
}

function min(value1,value2)
{
    return (value1 < value2 ? value1 : value2)
}

function abs(value)
{
  return (value<0?-value:value);
}

function distance(thisFilename,atomNumber1,atomNumber2) {
    x1=x[thisFilename,atomNumber1]
    y1=y[thisFilename,atomNumber1]
    z1=z[thisFilename,atomNumber1]
    x2=x[thisFilename,atomNumber2]
    y2=y[thisFilename,atomNumber2]
    z2=z[thisFilename,atomNumber2]
    temp=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2)
    temp=sqrt(temp)
    return temp
}

# returns the distance between the centroid of atoms a[1-6] and b[1-6]
function centroidDistance(thisFilename, a1, a2, a3, a4, a5, a6, b1, b2, b3, b4, b5, b6) {
    x1 = (x[thisFilename,a1]+x[thisFilename,a2]+x[thisFilename,a3]+x[thisFilename,a4]+x[thisFilename,a5]+x[thisFilename,a6])/6.0    
    y1 = (y[thisFilename,a1]+y[thisFilename,a2]+y[thisFilename,a3]+y[thisFilename,a4]+y[thisFilename,a5]+y[thisFilename,a6])/6.0    
    z1 = (z[thisFilename,a1]+z[thisFilename,a2]+z[thisFilename,a3]+z[thisFilename,a4]+z[thisFilename,a5]+z[thisFilename,a6])/6.0    
    x2 = (x[thisFilename,b1]+x[thisFilename,b2]+x[thisFilename,b3]+x[thisFilename,b4]+x[thisFilename,b5]+x[thisFilename,b6])/6.0    
    y2 = (y[thisFilename,b1]+y[thisFilename,b2]+y[thisFilename,b3]+y[thisFilename,b4]+y[thisFilename,b5]+y[thisFilename,b6])/6.0    
    z2 = (z[thisFilename,b1]+z[thisFilename,b2]+z[thisFilename,b3]+z[thisFilename,b4]+z[thisFilename,b5]+z[thisFilename,b6])/6.0    
    temp=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2)
    temp=sqrt(temp)
    return temp
}

function angle(thisFilename,atomNumber1,atomNumber2,atomNumber3) {
    value=((-distance(thisFilename,atomNumber1,atomNumber3)^2+distance(thisFilename,atomNumber1,atomNumber2)^2+distance(thisFilename,atomNumber2,atomNumber3)^2)/(2*distance(thisFilename,atomNumber1,atomNumber2)*distance(thisFilename,atomNumber2,atomNumber3)))
   return acos(value)
}

function asin(x) { return (180/3.141592)*atan2(x, sqrt(1-x*x)) }

function acos(x) { return (180/3.141592)*atan2(sqrt(1-x*x), x) }

function atan(x) { return (180/3.141592)*atan2(x,1) }

function dihedral(thisFilename,atomNumber1,atomNumber2,atomNumber3,atomNumber4) {
    x1=x[thisFilename,atomNumber1]
    y1=y[thisFilename,atomNumber1]
    z1=z[thisFilename,atomNumber1]
    x2=x[thisFilename,atomNumber2]
    y2=y[thisFilename,atomNumber2]
    z2=z[thisFilename,atomNumber2]
    x3=x[thisFilename,atomNumber3]
    y3=y[thisFilename,atomNumber3]
    z3=z[thisFilename,atomNumber3]
    x4=x[thisFilename,atomNumber4]
    y4=y[thisFilename,atomNumber4]
    z4=z[thisFilename,atomNumber4]
   
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

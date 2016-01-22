# takes the lowest gradient point and puts it into a file
# set filename and header in END block!

FNR == 1 {
    fileCount++
    filenames[fileCount]=FILENAME
    minRMSforce[fileCount]=-1
}

/Standard orientation/,/Rotational constants/ {
    if (match($0,"[a-zA-Z]") == 0 && NF==6)
        {
            thisIteration = iterations[fileCount]+1
            atomSymbol[$1]=$2
            x[fileCount,thisIteration,$1]=$4
            y[fileCount,thisIteration,$1]=$5
            z[fileCount,thisIteration,$1]=$6
        }
}

/SCF Done/ {
    if ( $5 < minEnergy || minEnergy == 0.0 )
        minEnergy = $5
    iterations[fileCount]++
    thisIteration = iterations[fileCount]
    energy[fileCount,thisIteration]=$5
}

/NAtoms=/ {
    if ( numberOfAtoms == 0 )
        {
            numberOfAtoms=$2
        }
    else if ( numberOfAtoms != $2 )
        {
            print "\nmismatch in number of atoms!", numberOfAtoms, $2, FILENAME
            exitCode=1
            exit
        }
}

/ Maximum Force/ && NF == 5 {
    thisIteration = iterations[fileCount]
    maxForces[fileCount,thisIteration]=$3
    getline
    RMSforces[fileCount,thisIteration]=$3
    if ( minRMSforce[fileCount] == -1 || $3 < minRMSforce[fileCount] )
        {
            minRMSforce[fileCount]=$3
            minRMSforceIteration[fileCount] = thisIteration
        }
    getline
    maxDisplacements[fileCount, thisIteration]=$3
    getline
    RMSdisplacements[fileCount, thisIteration]=$3
}

/Normal termination/ {
    finished[fileCount]++
}

/^ Frequencies --/ {
    didFrequencies[fileCount]=1
    for (i=3; i <= NF; i++)
        {
            if ($i < 0)
                {
                    if (fileCount in imaginaries)
                        imaginaries[fileCount]=imaginaries[fileCount] "," sprintf("%.0f",$i)
                    else
                        imaginaries[fileCount]=sprintf("%.0f",$i)
                }
        }
}

/Enter / && NF == 2 {
    numberOfFields=split($2,fields,"/")
    thisLink=fields[numberOfFields]
    gsub("[^0-9]","",thisLink)
    if (thisLink==9999)
        lastLink[fileCount]=""
    else if (thisLink==1002 || thisLink==1101 || thisLink==1110)
        lastLink[fileCount]="freq"
    else if (thisLink==502 || thisLink == 301)
        lastLink[fileCount]="scf"
    else if (thisLink==701 || thisLink==702 || thisLink==703)
        lastLink[fileCount]="grad"
    else
        lastLink[fileCount]=thisLink
}

END {
    printf "\n"
    if ( exitCode != 0 )
        exit
    
    outputCount = 0
    for (i=1; i <= fileCount; i++)
        {
            bestIteration = minRMSforceIteration[i]
            if ( bestIteration == 0 )
                continue
            outputCount++

            #outputFilename=sprintf("indocat_mannosyl_diax_anomer_BnOH_ts_%03d.gjf", outputCount)
            outputFilename=sprintf("indocat_mannosyl_diax_BnOH_ts_%03d.gjf", outputCount)
            print "%chk=checkpoint.chk\n%nprocshared=16\n%mem=32GB" > outputFilename
            printf "#p m062x/6-31g(d) scrf=(pcm,solvent=toluene) opt=(ts,calcfc,noeigentest,maxstep=5) freq=noraman pop=none\n\n" >> outputFilename
            printf "min force %.1f\n\n0 1\n", minRMSforce[i]*1E5 >> outputFilename
            printf "file %3d   %50s   min force %.1f\n", outputCount, filenames[i], minRMSforce[i]*1E5

            for (j=1; j < numberOfAtoms; j++)
                printf "%2s %12.8f %12.8f %12.8f\n", atomSymbol[j], x[i,bestIteration,j], y[i,bestIteration,j], z[i,bestIteration,j] >> outputFilename

            printf "\n\n" >> outputFilename
        }

    printf "%d files written\n", outputCount
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

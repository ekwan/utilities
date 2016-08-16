# takes a bunch of .out files and makes jobs where
# a bond distance is frozen, followed by a ts optimization
#
# if the job doesn't have "ts" in the title, a regular unconstrained
# optimization and frequency job will be run
#
# usage awk -f ts_theories.awk *.out

FNR == 1 {
    fileCount++
    filenames[fileCount]=FILENAME
}

/Standard orientation/,/Rotational constants/ {
    if ( NF == 6 && match($0,"[a-zA-Z]") == 0 )
        {
            symbol[fileCount,$1]=$2
            x[fileCount,$1]=$4
            y[fileCount,$1]=$5
            z[fileCount,$1]=$6
            atoms[fileCount]=$1
        }
}

END {
#    temp = "pbe0 pbe0_d3bj b3lyp b3lyp_d3bj b97d b97d3 bhandhlyp bhandh lcwpbe_d3bj bp86_d3bj bmk_d3bj cam_b3lyp_d3bj m062x m062x_d3 m06_hf_d3 m06l_hf_d3 hf"
#    temp = "blyp_d3bj tpsstpss_d3bj vsxc sogga11 b3p86 b98 b971 b972 o3lyp mpw1pw91 mpw1lyp mpw1pbe mpw3pbe m11 n12sx mn12sx b1b95 b3pw91 m06l_d3 b97d3_d3bj"
#    temp = "mp2"
#    temp = "ccsd"
    temp = "b3lyp m062x"
    split(temp, theory_abbreviations, " ")

    #temp = "pbe1pbe;pbe1pbe empiricaldispersion=gd3bj;b3lyp;b3lyp empiricaldispersion=gd3bj;b97d;b97d empiricaldispersion=gd3;bhandhlyp;bhandh;"
    #temp = temp "lc-wpbe empiricaldispersion=gd3bj;bp86 empiricaldispersion=gd3bj;bmk empiricaldispersion=gd3bj;cam-b3lyp empiricaldispersion=gd3bj;"
    #temp = temp "m062x;m062x empiricaldispersion=gd3;m06hf empiricaldispersion=gd3;m06l empiricaldispersion=gd3;hf"
    #temp = "blyp empiricaldispersion=gd3bj;tpsstpss empiricaldispersion=gd3bj;vsxc;sogga11;b3p86;b98;b971;b972;o3lyp;mpw1pw91;mpw1lyp;"
    #temp = temp "mpw1pbe;mpw3pbe;m11;n12sx;mn12sx;b1b95;b3pw91;m06l empiricaldispersion=gd3;b97d3 empiricaldispersion=gd3bj"
    #temp = "mp2"
    #temp = "ccsd"
    temp = "b3lyp;m062x"
    split(temp, theories, ";")

    if ( length(theory_abbreviations) != length(theories) )
        {
            print "level of theory sizes don't match"
            exit 1
        }

    for (i=1; i < length(theories); i++)
        printf "%30s %s\n", theory_abbreviations[i], theories[i]

    #basis_sets[1] = "6-31g*"
    #basis_sets[2] = "cc-pVDZ"
    #basis_sets[3] = "cc-pVTZ"
    basis_sets[4] = "cc-pVQZ"

    #basis_abbreviations[1] = "631gd"
    #basis_abbreviations[2] = "dz"
    #basis_abbreviations[3] = "tz"
    basis_abbreviations[4] = "qz"

    #solvents[1] = ""
    solvents[2] = "scrf=(solvent=o-xylene,pcm)"
    #solvents[3] = "scrf=(solvent=o-xylene,smd)"

    #solvent_abbreviations[1] = "gas"
    solvent_abbreviations[2] = "pcm"
    #solvent_abbreviations[3] = "smd"

    header = "%chk=checkpoint.chk\n%nprocshared=12\n%mem=12GB"
    vanilla_route   = "#p opt freq=noraman pop=none"
    standard_route1 = "#p opt=(modredundant,maxcycles=200) pop=none"
    standard_route2 = "#p opt=(ts,calcfc,noeigentest,maxstep=5,nofreeze,maxcycles=100) freq=noraman pop=none geom=allcheck guess=read"

    printf "" > "routes.txt"
    for (t in theories)
        {
            theory = theories[t]
            theory_abbreviation = theory_abbreviations[t]
            for (b in basis_sets)
                {
                    basis = basis_sets[b]
                    basis_abbreviation = basis_abbreviations[b]

                    for (s in solvents)
                        {
                            solvent = solvents[s]
                            solvent_abbreviation = solvent_abbreviations[s]
                            method = theory " " basis " " solvent
                            method_abbreviation = theory_abbreviation "-" basis_abbreviation "-" solvent_abbreviation
                            printf "#p %s\n", method >> "routes.txt"
                            for (i=1; i <= fileCount; i++)
                                {
                                    filename = filenames[i]
                                    gsub(".out","",filename)
                                    filename = "gjf/" filename "-" method_abbreviation ".gjf"
                                    print filename
                                    is_ts = match(filename, "endo") > 0 ? 1 : 0
                                    if ( is_ts == 0 )
                                        {
                                            print header "\n" vanilla_route " " method > filename
                                            print "\ntitle\n\n0 1" >> filename
                                        }
                                    else
                                        {
                                            print header "\n" standard_route1 " " method > filename
                                            print "\ntitle\n\n0 1" >> filename
                                        }

                                    for (j=1; j <= atoms[i]; j++)
                                        printf "%2d %15.8f %15.8f %15.8f\n", symbol[i,j], x[i,j], y[i,j], z[i,j]  >> filename
                                    if ( is_ts == 1 )
                                        print "\nB 2 6 F\nB 4 5 F\n\n--Link1--\n" header "\n" standard_route2 " " method >> filename
                                    printf "\n\n" >> filename
                                }
                        }
                }
        }
}

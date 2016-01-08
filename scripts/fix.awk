FNR > 1 {
    gsub("opt", "opt=modredundant", $0)
    gsub("midix", "midix scrf=(pcm,solvent=toluene)")
    print
}

/^ 167$/ {
    print "\nB 119 138 F\n\n"
}

# converts a multi-structure gjf to a multi-structure xyz
# assumes no repeated blank lines

{
    stripped = strip($0)
    if (length(stripped) == 0) {
        blanks++
    }

    if (blanks == 1) {
        title = title $0
    }
    else if (blanks == 2 && NF==4) {
        n++
        geometry[n]=$0
    }
    else if (blanks > 2) {
        print(n)
        print(title)
        for (i=1; i<=n; i++)
            print(geometry[i])
        n=0
        title=""
        blanks=0
    }
}

function ltrim(s) { sub(/^[ \t\r\n]+/, "", s); return s }
function rtrim(s) { sub(/[ \t\r\n]+$/, "", s); return s }
function strip(s) { return rtrim(ltrim(s)); }

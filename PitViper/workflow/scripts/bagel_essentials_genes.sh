

gawk '($2 > 0) {print $0}' $1 > $2

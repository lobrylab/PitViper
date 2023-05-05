# Description: Extracts essential genes from Bagel output

gawk '($2 > 0) {print $0}' $1 > $2

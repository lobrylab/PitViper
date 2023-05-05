# Convert a library file to a fasta file

awk -F ',' '{print ">"$1"\n"$2}' $1 > $2

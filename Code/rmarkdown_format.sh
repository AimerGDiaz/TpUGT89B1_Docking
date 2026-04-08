
awk  'NR == 1 {for (i = 1; i <= NF; i++) printf $i" | "; print ""; for (i = 1; i <= NF; i++) printf " --- |";print "";}NR >1 {gsub(/.*_/,"",$1);printf "model_"$1" | "; for (i = 2; i <= NF; i++) printf("%0.2f | ",$i); print ""}' $1

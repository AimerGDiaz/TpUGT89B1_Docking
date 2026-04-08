grep   "\:" *_summary_confidence*.json | awk -F':' '{gsub(/.*nces_|.json/,"",$1);print "model_"$1 , $2 , $3} '  | grep "\." | grep -v "\["  > test.txt
awk '
{
    gsub(/[,"]/, "");   # Remove quotes and commas
    model=$1; metric=$2; value=$3;
    data[model][metric] = value;  
    metrics[metric] = 1;  # Track all unique column headers
}
END {
    # Print header row
    printf "models";
    for (m in metrics) printf "\t%s", m;
    print "";  # New line

    # Print each model and its corresponding values
    for (mod in data) {
        printf "%s", mod;
        for (m in metrics) printf "\t%s", (data[mod][m] != "" ? data[mod][m] : "NA"); 
        print "";
    }
}' test.txt > models_ranking.txt 

rm test.txt 
cat models_ranking.txt


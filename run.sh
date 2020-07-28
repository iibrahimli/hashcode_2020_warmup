#! /bin/bash

total=$((0))
for file in $(ls in/); do
    score=$(./main "in/$file" "out/${file%.*}.out" | grep "score" | cut -d' ' -f2)
    printf "%-16s: %12s\n" "$file" "$score"
    total=$((total+score))
done
echo "------------------------------"
printf "%-16s: %12s\n" "total score" "$total"

echo "sum of cul 2"
awk '{sum += $2}; END {print sum}' tigr6.chrlen
echo "seprate by \t, default is space,\s"
awk -F"\t" '{sum += $2}; END {print sum}' tigr6.chrlen


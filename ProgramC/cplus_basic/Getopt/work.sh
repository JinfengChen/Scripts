echo "getopt with parameter parse by getopt"
g++ -o getopt getopt.cpp
./getopt -o chr01 -f 1 -m 200
./getopt

echo "no parse of parameter"
g++ -o argv argv.cpp
./argv -o chr01 -f 1 -m 200
./argv


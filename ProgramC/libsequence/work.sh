g++ -o revcom revcom.cc -lsequence
g++ -o valid valid.cc -lsequence
g++ -o alignment alignment.cc -lsequence

gcc -o valid_dna valid_dna.cc -lsequence -lboost_regex

gcc -o convert convert.cpp -lsequence



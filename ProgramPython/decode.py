import sys
import re

def main(code):
    assci = range(97,123)
    assci.extend([97,98])
    with open (code,'r') as codefh:
        for line in codefh:
            line.rstrip()
            words = line.split(' ')
            decode = []
            for w in words:
                newword = []
                for i in range(len(w)):
                    letter=w[i]
                    if (ord(letter) <= 122  and ord(letter) >= 97):
                        new = unichr(assci[assci.index(int(ord(letter)))+2])
                        newword.append(new)
                    else:
                        new = letter
                        newword.append(new)
                decode.append("".join(newword))
            print " ".join(decode)

if __name__ == '__main__':
    if len (sys.argv) < 2:
        print 'No parameter specified'
        sys.exit(2)
    main(sys.argv[1])




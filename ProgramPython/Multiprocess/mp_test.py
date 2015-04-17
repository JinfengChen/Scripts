import multiprocessing as mp
import random
import string

# Define an output queue
output = mp.Queue()

# define a example function
def rand_string(length, output):
    """ Generates a random string of numbers, lower- and uppercase chars. """
    rand_str = ''.join(random.choice(
                    string.ascii_lowercase
                    + string.ascii_uppercase
                    + string.digits)
               for i in range(length))
    ofile = open('%s.txt' %(rand_str), 'w')
    print >> ofile, rand_str
    ofile.close()
    output.put(rand_str)

# Setup a list of processes that we want to run
#processes = [mp.Process(target=rand_string, args=(5, output)) for x in range(4)]
processes = []
for x in range(4):
    p = mp.Process(target=rand_string, args=(5, output))
    processes.append(p)
    p.start()
 
# Run processes
#for p in processes:
#    p.start()

# Exit the completed processes
for p in processes:
    p.join()

# Get process results from the output queue
#results = [output.get() for p in processes]
#for p in processes:
#    print output.get()

#print(results)

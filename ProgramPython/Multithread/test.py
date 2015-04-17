import sys
sys.path.append("/rhome/cjinfeng/software/ProgramPython/module/lib")
from multi import worker

def main():
    #queue = Queue.Queue()
    cmd = ['python write.py 1',
    'python write.py 2',
    'python write.py 3',
    'python write.py 4',
    'python write.py 5',
    'python write.py 6',
    'python write.py 7']

    worker(cmd)
    '''Create a thread pool with queue'''
    #for i in range(5):
    #    t = Runner(queue)
    #    t.daemon = True
    #    t.start()
    
    ''''''
    #for c in cmd:
        #print c
    #    queue.put(c)
   
    #queue.join()
    print 'All job done'

'''
    runs = []
    for c in cmd:
        run = Worker(c)
        print 'Run one'
        run.start()
        runs.append(run)
    
    for r in runs:
        r.join()
'''

if __name__ == "__main__":
    main()


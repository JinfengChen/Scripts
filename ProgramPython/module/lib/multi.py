import os
import Queue
import threading
import time

class Runner(threading.Thread):
    '''Running cmd using multi-thread'''
    def __init__(self, queue):
        threading.Thread.__init__(self)
        self.queue = queue

    def run(self):
        '''run will be called when the object was created'''
        while True:
            cmd = self.queue.get()
            self.runner(cmd)
            self.queue.task_done()


    def runner(self, cmd):
        '''run cmdline'''
        os.system(cmd)

def worker(cmd, run=5):
    queue = Queue.Queue()
    for i in range(run):
        t = Runner(queue)
        t.daemon = True
        t.start()
    for c in cmd:
        queue.put(c)
    queue.join()

def main():
    '''
    if use outside this file, "from multi import worker"
    '''
    cmd = ['python write.py 1',
    'python write.py 2',
    'python write.py 3',
    'python write.py 4',
    'python write.py 5',
    'python write.py 6',
    'python write.py 7']

    worker(cmd)
    print 'All job done'

if __name__ == "__main__":
    main()


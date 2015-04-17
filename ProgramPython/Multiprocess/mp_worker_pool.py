#
# A test of `multiprocessing.Pool` class
#
# Copyright (c) 2006-2008, R Oudkerk
# All rights reserved.
#

import multiprocessing
import time
import random
import sys

#
# Functions used by test code
#

def calculate(func, args):
    result = func(*args)
    return '%s says that %s%s = %s' % (
        multiprocessing.current_process().name,
        func.__name__, args, result
        )

def calculatestar(args):
    return calculate(*args)

def mul(a, b):
    time.sleep(0.5*random.random())
    return a * b

def plus(a, b):
    time.sleep(0.5*random.random())
    return a + b

def f(x):
    return 1.0 / (x-5.0)

def pow3(x):
    return x**3

def noop(x):
    pass

#
# Test code
#

def test():
    print 'cpu_count() = %d\n' % multiprocessing.cpu_count()

    #
    # Create pool
    #

    PROCESSES = 4
    print 'Creating pool with %d processes\n' % PROCESSES
    pool = multiprocessing.Pool(PROCESSES)
    print 'pool = %s' % pool
    print

    #
    # Tests
    #

    TASKS = [(mul, (i, 7)) for i in range(10)] + \
            [(plus, (i, 8)) for i in range(10)]

    results = [pool.apply_async(calculate, t) for t in TASKS]
    imap_it = pool.imap(calculatestar, TASKS)
    imap_unordered_it = pool.imap_unordered(calculatestar, TASKS)

    print 'Ordered results using pool.apply_async():'
    for r in results:
        print '\t', r.get()
    print

    print 'Ordered results using pool.imap():'
    for x in imap_it:
        print '\t', x
    print

    print 'Unordered results using pool.imap_unordered():'
    for x in imap_unordered_it:
        print '\t', x
    print

    print 'Ordered results using pool.map() --- will block till complete:'
    for x in pool.map(calculatestar, TASKS):
        print '\t', x
    print

    #
    # Simple benchmarks
    #

    N = 100000
    print 'def pow3(x): return x**3'

    t = time.time()
    A = map(pow3, xrange(N))
    print '\tmap(pow3, xrange(%d)):\n\t\t%s seconds' % \
          (N, time.time() - t)

    t = time.time()
    B = pool.map(pow3, xrange(N))
    print '\tpool.map(pow3, xrange(%d)):\n\t\t%s seconds' % \
          (N, time.time() - t)

    t = time.time()
    C = list(pool.imap(pow3, xrange(N), chunksize=N//8))
    print '\tlist(pool.imap(pow3, xrange(%d), chunksize=%d)):\n\t\t%s' \
          ' seconds' % (N, N//8, time.time() - t)

    assert A == B == C, (len(A), len(B), len(C))
    print

    L = [None] * 1000000
    print 'def noop(x): pass'
    print 'L = [None] * 1000000'

    t = time.time()
    A = map(noop, L)
    print '\tmap(noop, L):\n\t\t%s seconds' % \
          (time.time() - t)

    t = time.time()
    B = pool.map(noop, L)
    print '\tpool.map(noop, L):\n\t\t%s seconds' % \
          (time.time() - t)

    t = time.time()
    C = list(pool.imap(noop, L, chunksize=len(L)//8))
    print '\tlist(pool.imap(noop, L, chunksize=%d)):\n\t\t%s seconds' % \
          (len(L)//8, time.time() - t)

    assert A == B == C, (len(A), len(B), len(C))
    print

    del A, B, C, L

if __name__ == '__main__':
    multiprocessing.freeze_support()
    test()

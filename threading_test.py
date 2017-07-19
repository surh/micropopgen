#!/usr/bin/env python3
from functools import partial
from itertools import repeat
from multiprocessing import Pool, freeze_support

def func(a, b):
    return a['1'] + b

def main():
    a_args = [{'1' : 1},{'1' : 2},{'1' : 3}]
    second_arg = 1
    with Pool(2) as pool:
        L = pool.starmap(func, [({'1' : 1}, 1), ({'1' : 2}, 1), ({'1' : 3}, 1)])
        M = pool.starmap(func, zip(a_args, repeat(second_arg)))
        N = pool.map(partial(func, b=second_arg), a_args)
        print(L)
        print(M)
        print(N)
        assert L == M == N

if __name__=="__main__":
    #freeze_support()
    main()


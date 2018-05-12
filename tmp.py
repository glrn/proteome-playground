from multiprocessing import Pool
from time import sleep

def f(args):
    a = args[0]
    b = args[1]
    c = args[2]
    return a+b+c

if __name__ == '__main__':
    p = Pool(2)
    print(p.map(f, [[1,2,3],[4,5,6]]))
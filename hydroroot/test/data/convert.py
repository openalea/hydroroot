import pandas as pd

def convert_line(line):
    return '\t'.join(line.strip().split())


def convert(fn1):
    f = open(fn1)
    lines = [convert_line(l) for l in f.readlines()]
    content = '\n'.join(lines)+'\n'
    f.close()

    f2 = open(fn1,'w')
    f2.write(content)
    f2.close()

    return True



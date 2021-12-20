#!/usr/bin/python


import argparse
import string
import os


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('--data', required=True, help="You MUST provide the name of data file for filter step.")
    parser.add_argument('--threshold', default="20", help="The mapping reads number threshold under which "
                                                          "those genes will be filtered out.")
    parser.add_argument('-d', '--dir', default=".", help="The path of output directory.")
    parser.add_argument('--prefix', default='', help="The prefix for all the output files.")
    args = parser.parse_args()

    count = 1
    repository = []
    scale = []
    whether_first_line = 1
    pattern = "JUST.FOR.START"

    f = open(args.data, 'r')
    w = open(os.join.path(args.dir, 'data_filter'), 'w')
    c = open(os.join.path(args.dir, 'data_count'), 'w')

    line = f.readline()
    while line:
        if whether_first_line:
            w.write(line)
            whether_first_line = 0
            line = f.readline()
            continue

        com = line.split(' ')

        if pattern != com[1] and count == 1:
            pattern = com[1]
            repository = [line]
        elif pattern == com[1]:
            repository.append(line)
            count += 1
        else:
            for i in range(0, count):
                tmp = repository[i].split(' ')
                print(tmp)
                tmp = [string.atoi(x) for x in tmp if x.isdigit()]

                scale.append(max(tmp[1:]))
            if sum(scale) > max(scale) > args.threshold:
                for i in repository:
                    w.write(i)
                c.write(str(count) + "\n")
            count = 1
            pattern = com[1]
            repository = [line]
            scale = []
        line = f.readline()

    f.close()
    w.close()
    c.close()


if __name__ == '__main__':
    main()
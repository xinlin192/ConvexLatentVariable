
import csv

with open('wholesale.csv') as csvfile:
    wholesale_reader = csv.reader(csvfile, delimiter=',')
    line_num = 0
    dtable = {}
    dmax = 0
    for line in wholesale_reader:
        if line_num == 0:
            line_num += 1
            continue
        d = (line[0], line[1])
        if not dtable.has_key(d):
            dmax += 1
            dindex = dmax
            dtable.update({d:dindex})
        dindex = dtable[d]
        outstr = ""
        outstr += str(dindex)
        for idx in range(2, 8):
            outstr += " " + str(idx-1)+":" + line[idx]
        print outstr
        line_num += 1

import os

input = open("test.faa")

tcdbSystems = {}

for line in input:
    if(">" in line):
        first = line.split("-")[0][1:]
        if(first in tcdbSystems.keys()):
            tcdbSystems.get(first).append(line.strip(">\n"))
        else:
            tcdbSystems[first] = [line.strip(">\n")]

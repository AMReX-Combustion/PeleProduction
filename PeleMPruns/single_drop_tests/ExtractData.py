
import os
import csv
import numpy as np
from CaseInfo import *

def ExtractData(case, outfile):
    pltfiles = []
    sprayfiles = []
    for x in os.listdir("./" + case.name + "/"):
        if (x.startswith("plt")):
            pltfiles.append(case.name + "/" + x)
        if (x.endswith(".p3d")):
            sprayfiles.append(case.name + "/" + x)
    def get_step(fn):
        x = fn.split()
        res = []
        for i in x:
            if (i.isnumeric()):
                res.append(i)
        return [fn, res]
    pltfiles = sorted(pltfiles, key=get_step)
    sprayfiles = sorted(sprayfiles, key=get_step)
    alltime = []
    # Open Header in the plt files and extract the solution time
    for cf in pltfiles:
        curfile = cf + "/Header"
        Lines = []
        with open(curfile, 'r') as fn:
            Lines = fn.readlines()
            numcomp = int(Lines[1])
            timeline = numcomp + 3
        alltime.append(float(Lines[timeline]))
    # Column designations in the spray*.p3d files
    numspec = 2 # Liquid fuel components
    dims = 2 # Solution dimensions
    loccols = dims - 1
    velcols = loccols + 1
    tcol = velcols + dims
    dcol = tcol + 1
    mfcol = dcol + 1
    numcol = mfcol + numspec
    timevals = []
    vals = np.zeros([len(sprayfiles), numcol])
    crow = 5
    for i, cf in enumerate(sprayfiles):
        with open(cf, "r") as fn:
            Lines = fn.readlines()
            if (len(Lines) >= crow+1):
                timevals.append(alltime[i])
                sline = Lines[crow].split()
                for col in range(numcol):
                    vals[i][col] = float(sline[col])
    xconv = case.xconv
    yconv = case.yconv
    yexp = case.yexp
    modvals = []
    with open(outfile, "w+") as new_file:
        new_file.write("t, T, dd0, Y1, Y2\n")
        csv_writer = csv.writer(new_file, delimiter=',', lineterminator='\n')
        for k, tv in enumerate(timevals):
            dia = vals[k][dcol]
            dd0 = (dia * yconv)**yexp
            T = vals[k][tcol]
            Y1 = vals[k][mfcol]
            Y2 = vals[k][mfcol+1]
            outvals = [tv*xconv, T, dd0, Y1, Y2]
            modvals.append(outvals)
            csv_writer.writerow(outvals)
    modvals = np.array(modvals)
    return modvals

def ExtractRefVals(case):
    ldir = "ref_files/" + case.name + "/"
    fnames = ["refdvals.csv", "refTvals.csv", "refYvals.csv"]
    reffiles = []
    def getdata(fname):
        vals = []
        with open(fname, 'r') as rf:
            line0 = rf.readline()
            lines = rf.readlines()
            for line in lines:
                sline = line.split(',')
                ovals = [float(sline[0]), float(sline[1])]
                vals.append(ovals)
        vals = np.array(vals)
        return vals
    dvals = None
    tvals = None
    yvals = None
    cname = fnames[0]
    if (cname in os.listdir(ldir)):
        dvals = getdata(ldir + cname)
    cname = fnames[1]
    if (cname in os.listdir(ldir)):
        tvals = getdata(ldir + cname)
    cname = fnames[2]
    if (cname in os.listdir(ldir)):
        yvals = getdata(ldir + cname)
    return [dvals, tvals, yvals]

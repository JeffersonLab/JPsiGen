#!/usr/bin/env python2
#encoding: UTF-8

import getopt
import os
import sys

def usage():
    print("The usage:");
    print("python GenerateGenOptions.py -n Nsim -e Ebeam -t t_lim --Egmin Egmin --Egmax Egmax --q2Cut Q2CutValue -l (write output in LUND file)");
    print("As an example ");
    print("python GenerateGenOptions.py -n 100000 --LUND --Nperfile 35000 -e 10.6 -t -6 --q2Cut 0.02  --ltarg 5 --Run");

def main():

# Default values of generator parameters
    Nsim = 100000
    NPerFile = 10000
    Eb = 10.6
    tLim = -6.
    Ltarg = 5.
    EgMin = 4
    EgMax = 10.6
    Q2Cut = 0.02
    LUND = 0
    Run = False;

    try:
        #opts, args = getopt.getopt(sys.argv[1:], "n:e:t:lhor", ["help", "output=", "Nperfile=", "Egmin=", "Egmax=", "ltarg=", "q2Cut=", "LUND", "Run"])
        opts, args = getopt.getopt(sys.argv[1:], "n:e:t:lhor", ["help", "output=", "Nperfile=", "Egmin=", "Egmax=", "ltarg=", "q2Cut=", "LUND", "Run"])
    except getopt.GetoptError as err:
        # print help information and exit:
        print str(err)  # will print something like "option -e not recognized"
        usage()
        sys.exit(2)
    for o, a in opts:
        if o == "-n":
            Nsim = a
        elif o in("--Nperfile"):
            NPerFile = a
        elif o == "-e":
            Eb = a
        elif o == "-t":
            if(float(a) > 0):
                print("You provided tMax as " + a)
                print("tMax should be a negative number.")
                print("Exiting")
                sys.exit()               
            tLim = a
        elif o in ("--Egmin"):
            EgMin = a;
        elif o in ("--Egmax"):
            EgMax = a;
        elif o in ("--ltarg"):
            Ltarg = a;
        elif o in ("--q2Cut"):
            Q2Cut = a;
        elif o in ("-l", "--LUND"):
            LUND = 1
            print "KUKU"
            print "LUND is " + str(LUND)
        elif o in ("-r", "--Run"):
            Run = True;
        elif o in ("-h", "--help"):
            usage()
            sys.exit()
        else:
            assert False, "unhandled option"


    optFile = open("GenOptions.dat", "w")
    
    optFile.write("Nsim     " + str(Nsim) + "\n");
    optFile.write("NPerFile " + str(NPerFile) + "\n");
    optFile.write("Eb       " + str(Eb) + "\n");
    optFile.write("tLim     " + str(tLim)  + "\n");
    optFile.write("lTarg    " + str(Ltarg) + "\n");
    optFile.write("EgMin    " + str(EgMin) + "\n");
    optFile.write("EgMax    " + str(EgMax) + "\n");
    optFile.write("Q2Cut    " + str(Q2Cut) + "\n");
    optFile.write("LUND     " + str(LUND));

    optFile.close()
    
    
    if(Run):
        Cmd = "./JPsiGen.exe"
        print("Running the generator from the python code \n")
        os.system(Cmd);

if __name__ == "__main__":
    main();
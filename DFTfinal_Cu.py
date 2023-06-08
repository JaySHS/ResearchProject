from pymatgen.io import pwscf
import pymatgen.core as core
import numpy as np
import subprocess
import matplotlib.pyplot as plt
from pathlib import Path
import re
import csv
import time
import os
import eoslib as eos

def Relaxcheck(v, r):
    v_list = []
    if(os.path.isfile(alloyfolder+"C16_data.csv")):
        with open(alloyfolder+"C16_data.csv") as f:
            c16line=f.readlines()
            for line in c16line:
                v_list.append(float(line.split(',')[0]))
            f.close()
        for i in range(len(v_list)):
            if(v_list[i]==v):
                return False
        if(os.path.isfile(infolder+"C16_v%.2f"%v+"_ca0.70_relax"+str(r+1)+".in")):
            return False
        else:
            return True
    else:
        return True
    

def debug(cmd):
    with open(alloyfolder+"progress.dat", 'a') as f:
        f.write(cmd+"\n")
        f.close()

def allscfdone(path, v, rindex):
    with open(path+'C16v%.2f_relax%.2fscfdone.dat'%(v, rindex), 'r') as f:
        check = f.read()
        return (check == 'xxxxxxxx')

def checkoccupy(path):
    while(os.path.isfile(path+"Occupied")):
        time.sleep(1)

def occupy(path):
    with open(path+"Occupied", 'x') as f:
        f.close()

def deoccupy(path):
    os.remove(path+"Occupied")

def bash_run(cmd):
    subprocess.run(cmd, shell=True, executable='/bin/bash')

def checkdir(dirname):
    try:
        subprocess.run(["mkdir", dirname])
    except:
        pass

def finalinput_c1(v):
    a = (v*4)**(1/3)
    structure = core.Structure(lattice=[[a, a, 0],[0, a, a], [a, 0, a]],
                                 species=[element1, element2, element2], 
                                 coords=[[0, 0, 0], [0.25, 0.25, 0.25], [-0.25, -0.25, -0.25]])
    pwinput = pwscf.PWInput(control = {"calculation": "scf", "pseudo_dir":"./pseudos/"},
                              system= {"ibrav":0, "ecutwfc":cut, "nbnd":40, "nat":3, 
                                       "occupations":"smearing", "smearing":"gaussian", "degauss":0.01},
                              electrons={"mixing_beta":0.6},
                              structure=structure,
                              pseudo = {element1: pseudodict[element1], element2: pseudodict[element2]},
                              kpoints_grid=(k, k, k), kpoints_shift=(1, 1, 1))
    return str(pwinput)

def finalinput_c16(v, caratio):
    a = (v*2/caratio)**(1/3)
    c = a*caratio
    structure = core.Structure(lattice=[[a/2, -a/2, c/2],[a/2, a/2, c/2], [-a/2, -a/2, c/2]],
                               species=[element1, element1, element2, element2, element2, element2], 
                               coords=relaxcoord)
    pwinput = pwscf.PWInput(control = {"calculation": "scf", "pseudo_dir":"./pseudos/"},
                            system= {"ibrav":0, "ecutwfc":cut, "nbnd":40, "nat":6},
                            electrons={"mixing_beta":0.6},
                            structure=structure,
                            pseudo = {element1: pseudodict[element1], element2: pseudodict[element2]},
                            kpoints_grid=(k, k, int(k/caratio)))
    return str(pwinput)

def odd(x):
    ix = int(x)
    if(ix%2 == 0):
        return ix
    else:
        return (ix-1)

pseudodict = {"Ag":"Ag.pbe-n-rrkjus_psl.1.0.0.UPF",
              "Au":"Au.pbe-n-rrkjus_psl.1.0.0.UPF",
              "Cu":"Cu.pbe-dn-rrkjus_psl.1.0.0.UPF",
              "In":"In.pbe-dn-rrkjus_psl.1.0.0.UPF"}

defcoord = [[0, 0.25, 0.25], [0, -0.25, -0.25], 
            [-0.5, 0.6304, -0.1304], [-0.5, 0.3696, 0.1304], 
            [0.239201, 0.1304, -0.3696], 
            [0.760799, -0.1304, -0.6304]]

cut = 100
k = 10
carange = np.linspace(0.7, 0.9, 8)
v0 = 700.0
dv = 20.0
vconvlim = 0.5
vconv = False

sep = 'CELL_'

folder = "./DFTData/"

element1 = "Cu"
element2 = "In"
alloyname = element1+element2+"2"

alloyfolder = folder+alloyname+"/"
infolder = alloyfolder+"in/"
outfolder = alloyfolder+"out/"
relaxfolder = alloyfolder+"relax/"

checkdir(folder)
checkdir(alloyfolder)
checkdir(infolder)
checkdir(outfolder)
checkdir(relaxfolder)

c1_v_data = []
c1_e_data = []

c16_v_data = []
c16_e_data = []

v = v0/2
while(not(vconv)):
    inputname = infolder+"C1_v%.2f.in"%v
    outputname = outfolder+"C1_v%.2f.out"%v
    if(not(os.path.isfile(inputname))):
        a = (v*4)**(1/3)
        structure = core.Structure(lattice=[[a, a, 0],[0, a, a], [a, 0, a]],
                                     species=[element1, element2, element2], 
                                     coords=[[0, 0, 0], [0.25, 0.25, 0.25], [-0.25, -0.25, -0.25]])
        pwinput = pwscf.PWInput(control = {"calculation": "scf", "prefix": alloyname+"_C1", "outdir":"./tmp/", 
                                             "pseudo_dir":"./pseudos/"},
                                  system= {"ibrav":2, "celldm(1)": a, "ecutwfc":cut, "nbnd":40, "nat":3, 
                                           "occupations":"smearing", "smearing":"gaussian", "degauss":0.01},
                                  electrons={"mixing_beta":0.6},
                                  structure=structure,
                                  pseudo = {element1: pseudodict[element1], element2: pseudodict[element2]},
                                  kpoints_grid=(k, k, k), kpoints_shift=(1, 1, 1))
        strin = str(pwinput)
        pwin, sepr, rm = strin.partition(sep)
    
        with open(inputname, "w+") as f:
            f.write(pwin)
            f.close()
        bash_run("mpirun -np 20 pw.x -i " + inputname + " |& tee " + outputname)
        debug("C1 calculation finished : V = %.2f"%v)

        output = pwscf.PWOutput(outputname)

        with open(alloyfolder+'C1_data.csv', 'a') as csvfile:
            csvwriter = csv.writer(csvfile)
            csvwriter.writerow([v, output.final_energy])
            csvfile.close()

    if(os.path.isfile(alloyfolder+'C1_data.csv')):
        checkoccupy(alloyfolder)
        occupy(alloyfolder)
        with open(alloyfolder+'C1_data.csv', 'r') as c1csv:
            c1line = c1csv.readlines()
            for line in c1line:
                c1data = line.split()[0].split(',')
                c1_v_data.append(float(c1data[0]))
                c1_e_data.append(float(c1data[1]))
            c1csv.close()
        deoccupy(alloyfolder)
    try:
        vol_c1, bmod_c1, p_c1, res_c1 = eos.eos_fit('birch', c1_v_data, c1_e_data, units=('Ryd', 'bohr'))
        ene_c1 = p_c1[0]
        
        if(np.abs(vol_c1-v)<vconvlim):
            vconv = True
        else:
            v = vol_c1

        checkoccupy(alloyfolder)
        occupy(alloyfolder)
        with open(alloyfolder+'c1v0.dat', 'w+') as f:
            f.write(str(vol_c1))
            f.close()
        deoccupy(alloyfolder)
    
    except:
        if(os.path.isfile(alloyfolder+'c1v0.dat')):
            checkoccupy(alloyfolder)
            occupy(alloyfolder)
            with open(alloyfolder+'c1v0.dat', 'r') as f:
                v_0 = float(f.read())
                f.close()
            if(np.abs(v_0-v)<vconvlim):
                vconv = True
            else:
                v = v_0
            deoccupy(alloyfolder)
        else:
            v = v+dv
    c1_v_data.clear()
    c1_e_data.clear()
    

input = finalinput_c1(vol_c1)
with open(alloyfolder+"C1final.in", "w+") as f:
    f.write(input)
    f.close()

debug("C1 calculation finished.")
debug("C1 minimum volume : %.3f"%vol_c1)
debug("C1 minimum energy : %.3f"%ene_c1)
debug('')

vconv = False
v = v0

emintmp = 0
ratiotmp = 0
while(not(vconv)):
    rindex = 0
    relaxcoord = defcoord
    relax = False
    while(not(relax)):
        emin = 0
        min_caratio = 0
        for caratio in carange:
            inputname = infolder+"C16_v%.2f"%v+"_ca%.2f_relax"%caratio+str(rindex)+".in"
            outputname = outfolder+"C16_v%.2f"%v+"_ca%.2f_relax"%caratio+str(rindex)+".out"
            try:
                with open(inputname) as f:
                    f.close()
            except:
                a = (v*2/caratio)**(1/3)
                c = a*caratio
                structure = core.Structure(lattice=[[a/2, -a/2, c/2],[a/2, a/2, c/2], [-a/2, -a/2, c/2]],
                                           species=[element1, element1, element2, element2, element2, element2], 
                                           coords=relaxcoord)
                pwinput = pwscf.PWInput(control = {"calculation": "scf", "prefix": alloyname+"_C16", "outdir":"./tmp/", 
                                                   "pseudo_dir":"./pseudos/"},
                                        system= {"ibrav":7, "celldm(1)": a, "celldm(3)": caratio, "ecutwfc":cut, "nbnd":40, "nat":6},
                                        electrons={"mixing_beta":0.6},
                                        structure=structure,
                                        pseudo = {element1: pseudodict[element1], element2: pseudodict[element2]},
                                        kpoints_grid=(k, k, k))
                strin = str(pwinput)
                pwin, sepr, rm = strin.partition(sep)
            
                with open(inputname, "w+") as f:
                    f.write(pwin)
                    f.close()
                bash_run("mpirun -np 20 pw.x -i " + inputname + " |& tee " + outputname)
                output = pwscf.PWOutput(outputname)

                emindat = alloyfolder+'C16eminv%.2f_relax%.2f.csv'%(v, rindex)

                if(os.path.isfile(emindat)):
                    checkoccupy(alloyfolder)
                    occupy(alloyfolder)
                    with open(emindat) as f:
                        eminrow = f.readline().split(',')
                        emin = float(eminrow[0])
                        min_caratio = float(eminrow[1])
                        if(emin > output.final_energy):
                            emin = output.final_energy
                            min_caratio = caratio
                        f.close()
                    deoccupy(alloyfolder)
                else:
                    emin = output.final_energy
                    min_caratio = caratio
                
                with open(emindat, 'w+') as csvfile:
                    csvwriter = csv.writer(csvfile)
                    csvwriter.writerow([emin, min_caratio])
                
                with open(alloyfolder+'C16v%.2f_relax%.2fscfdone.dat'%(v, rindex), 'a') as f:
                    f.write('x')
                    f.close()

        emindat = alloyfolder+'C16eminv%.2f_relax%.2f.csv'%(v, rindex)
        if(emin == 0):
            checkoccupy(alloyfolder)
            occupy(alloyfolder)
            with open(emindat) as f:
                eminrow = f.readline().split(',')
                emin = float(eminrow[0])
                min_caratio = float(eminrow[1])
                f.close()
            deoccupy(alloyfolder)

        relaxinname = relaxfolder+"Relax_V%.2f_%i.in"%(v, rindex)
        relaxoutname = relaxfolder+"Relax_V%.2f_%i.out"%(v, rindex)
        if(os.path.isfile(relaxinname)):
            while(Relaxcheck(v, rindex)):
                time.sleep(2)
        else:
            while(not(allscfdone(alloyfolder, v, rindex))):
                time.sleep(5)
            relax_a = (v*2/min_caratio)**(1/3)
            relax_c = relax_a*(min_caratio)

            rstructure = core.Structure(lattice=[[relax_a/2, -relax_a/2, relax_c/2],[relax_a/2, relax_a/2, relax_c/2],
                                                 [-relax_a/2, -relax_a/2, relax_c/2]],
                                            species=[element1, element1, element2, element2, element2, element2], 
                                            coords=relaxcoord)
            relaxinput = pwscf.PWInput(control = {"calculation": "relax", "prefix": alloyname+"_Relax", "outdir":"./tmp/", 
                                                  "pseudo_dir":"./pseudos/", "etot_conv_thr":0.000001, "forc_conv_thr":0.01},
                                       system= {"ibrav":7, "celldm(1)": relax_a, "celldm(3)": min_caratio, "ecutwfc":cut, "nat":6},
                                       electrons={"conv_thr":0.000001},
                                       structure=rstructure,
                                       pseudo = {element1: pseudodict[element1], element2: pseudodict[element2]},
                                       kpoints_grid=(k, k, k))
            strin = str(relaxinput)
            pwin, sepr, rm = strin.partition(sep)
            pwin = re.sub("1d", "1e", pwin)
        
            with open(relaxinname, "w+") as f:
                f.write(pwin)
                f.close()
            bash_run("mpirun -np 20 pw.x -i " + relaxinname + " |& tee " + relaxoutname)
            with open(alloyfolder+"Relaxed", 'x') as f:
                f.close()

        checkoccupy(alloyfolder)
        occupy(alloyfolder)

        of = open(relaxoutname, "r")
        ostring = of.read()
        of.close()
        of = open(relaxoutname, "r") 

        bfgregex = "bfgs\sconverged\sin\s+\d+\sscf\scycles\sand\s+(\d+)\sbfgs\ssteps"
        if re.search(bfgregex, ostring) != None:
            if int(re.search(bfgregex, ostring).group(1)) != 0:
                relaxed = False

                relaxcoord.clear()

                posregex = "Begin\sfinal\scoordinates"
                foundpos = False

                while not(foundpos):
                    line = of.readline()
                    if re.search(posregex, line) != None:
                        foundpos = True
        
                of.readline()
                of.readline()
                for x in range(0, 6):
                    line = of.readline()
                    coordlist = re.findall("-*\d+.\d+", line)
                    for i in range(0,3):
                        coordlist[i] = float(coordlist[i])
                    relaxcoord.append(coordlist)
                debug("Relax unfinished.")
            else:
                relax = True
                debug("Relax finished.")
        else:
            debug("Relax error.")
            relax = True
        of.close()
        rindex += 1
        if(os.path.isfile(alloyfolder+"Relaxed")):
            os.remove(alloyfolder+"Relaxed")
        deoccupy(alloyfolder)
    
    try:
        output = pwscf.PWOutput(relaxoutname)
        e_final = output.final_energy
        truerelax = True
    except:
        e_final = emin
        truerelax = False

    if(truerelax):
        with open(alloyfolder+'C16_data.csv', 'a') as csvfile:
            csvwriter = csv.writer(csvfile)
            csvwriter.writerow([v, e_final])
            csvfile.close()
    
        debug("C16 calculation finished : V = %.2f"%v)

    if(os.path.isfile(alloyfolder+'C16_data.csv')):
        checkoccupy(alloyfolder)
        occupy(alloyfolder)
        with open(alloyfolder+'C16_data.csv', 'r') as c16csv:
            c16line = c16csv.readlines()
            for line in c16line:
                c16data = line.split()[0].split(',')
                c16_v_data.append(float(c16data[0]))
                c16_e_data.append(float(c16data[1]))
            c16csv.close()
        deoccupy(alloyfolder)
    try:
        vol_c16, bmod_c16, p_c16, res_c16 = eos.eos_fit('birch', c16_v_data, c16_e_data, units=('Ryd', 'bohr'))
        ene_c16 = p_c16[0]

        if(np.abs(vol_c16-(v))<vconvlim):
            vconv = True
        else:
            v = vol_c16
        checkoccupy(alloyfolder)
        occupy(alloyfolder)
        with open(alloyfolder+'c16v0.dat', 'w+') as f:
            f.write(str(vol_c16))
            f.close()
        deoccupy(alloyfolder)

    except:
        if(os.path.isfile(alloyfolder+'c16v0.dat')):
            checkoccupy(alloyfolder)
            occupy(alloyfolder)
            with open(alloyfolder+'c16v0.dat', 'r') as f:
                v_0 = float(f.read())
                if(np.abs(v_0-v)<vconvlim):
                    vconv = True
                else:
                    v = v_0
                f.close()
            deoccupy(alloyfolder)
        else:
            v = v+dv
    
    c16_v_data.clear()
    c16_e_data.clear()

input = finalinput_c16(v, min_caratio)
with open(alloyfolder+"C16final.in", "w+") as f:
    f.write(input)
    f.close()

debug("C16 calculation finished.")
debug("C16 minimum volume : %.3f"%vol_c16)
debug("C16 minimum energy : %.3f"%ene_c16)
debug('')

debug("Run complete.")
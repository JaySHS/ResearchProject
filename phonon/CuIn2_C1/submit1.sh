#!/bin/bash

#$ -N CuIn2_C1_phonon1
#$ -S /bin/bash
#$ -cwd
#$ -o CuC1-001.out
#$ -e job_output1.err
#$ -m e
#$ -M shinhon@oregonstate.edu
#$ -pe orte 20

pw.x -i CuC1-001.in
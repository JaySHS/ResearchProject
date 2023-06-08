#!/bin/bash

#$ -N CuIn2_C1_phonon2
#$ -S /bin/bash
#$ -cwd
#$ -o CuC1-002.out
#$ -e job_output2.err
#$ -m e
#$ -M shinhon@oregonstate.edu
#$ -pe orte 20

pw.x -i CuC1-002.in
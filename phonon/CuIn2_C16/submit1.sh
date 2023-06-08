#!/bin/bash

#$ -N CuIn2_C16_phonon1
#$ -S /bin/bash
#$ -cwd
#$ -o CuC16-001.out
#$ -e job_output1.err
#$ -m e
#$ -M shinhon@oregonstate.edu
#$ -pe orte 40

pw.x -i CuC16-001.in
#!/bin/bash

#$ -N CuIn2_C16_phonon3
#$ -S /bin/bash
#$ -cwd
#$ -o CuC16-003.out
#$ -e job_output3.err
#$ -m e
#$ -M shinhon@oregonstate.edu
#$ -pe orte 40

pw.x -i CuC16-003.in
#!/bin/bash

#$ -N CuIn2_C16_phonon2
#$ -S /bin/bash
#$ -cwd
#$ -o CuC16-002.out
#$ -e job_output2.err
#$ -m e
#$ -M shinhon@oregonstate.edu
#$ -pe orte 40

pw.x -i CuC16-002.in
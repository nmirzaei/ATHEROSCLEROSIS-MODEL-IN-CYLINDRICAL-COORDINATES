#!/bin/bash
#SBATCH --job-name=Cylinder_Isotropic_Growth
#SBATCH --output=my_job_op%j.txt
#SBATCH --partition=mathsci
#SBATCH --mem=100G
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=nmirzaei@udel.edu
#SBATCH --mail-type=END


vpkg_require singularity
vpkg_require fenics
Sexec python3 main_script.py

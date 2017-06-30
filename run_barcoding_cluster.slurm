#!/usr/bin/env bash

# This you first segment you only need to do once.

setpkgs -a cmake
wget https://github.com/RuleWorld/bionetgen/archive/BioNetGen-2.2.6-stable.tar.gz
tar xzf BioNetGen-2.2.6-stable.tar.gz
mv bionetgen-BioNetGen-2.2.6-stable/ bionetgen
cd $HOME/bionetgen/bng2/Network3
cp Makefile.cmake Makefile
make Makefile
cd $HOME
export BNGPATH=$HOME/bionetgen
mkdir git
cd git
git clone -b bng_simulator https://github.com/LoLab-VU/pysb.git

#### For the run scripts
# all the slurm stuff up here

#SBATCH --mail-user=<corey.hayford@vanderbilt.edu>
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=00:10:00
#SBATCH --mem=1G
#SBATCH --output=barcoding_testjob.out

module load Anaconda2/4.3.1
export BNGPATH="$HOME/BioNetGen-2.3.0"
export PYTHONPATH="$HOME/git/pysb"

cd /scratch/hayforc/
python barcoding_model_final.py
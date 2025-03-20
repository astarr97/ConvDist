#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH -p hbfraser
#SBATCH --mem=32GB

#Added so it isn't deleted
wget https://cgl.gi.ucsc.edu/data/cactus/447-mammalian-2022v1.hal
mv 447-mammalian-2022v1.hal hg38.447way.hal

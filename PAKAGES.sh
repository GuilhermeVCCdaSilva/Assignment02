conda update conda

# Install sra-tools from conda and create a env
source ~/miniconda3/etc/profile.d/conda.sh
conda create -n sra
conda activate sra
conda install sra-tools
conda deactivate

# Install ipyrad and create env
source ~/miniconda3/etc/profile.d/conda.sh
conda config --add channels conda-forge
conda create -n ipyrad
conda activate ipyrad

conda install ipyrad -c bioconda
conda install -c bioconda muscle=3.8.1551
conda install -c bioconda vsearch=2.19
conda deactivate

# Install structure_threader  and create env
source ~/miniconda3/etc/profile.d/conda.sh
conda config --add channels conda-forge
conda create -n structure python=3.7 r r-rcolorbrewer bioconductor-snprelate bioconductor-lfa r-devtools -c bioconda 
conda activate structure  
pip install structure_threader
conda deactivate 

# Install R from apt
sudo apt install r-base
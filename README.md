# Amino Acid Propensity

### For mmseq2

conda install -c conda-forge -c bioconda mmseqs2

### For dssp

git clone https://github.com/PDB-REDO/dssp.git
cd dssp
cmake -S . -B build -DBUILD_PYTHON_MODULE=ON
cmake --build build
sudo cmake --install build

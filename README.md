# Protein-Internal-Energy
Final project for upper division course ECS 129   
By Kaustubh Deshpande  
https://web.cs.ucdavis.edu/~koehl/Teaching/ECS129/projects.html  
https://web.cs.ucdavis.edu/~koehl/Teaching/ECS129/Projects/Assignment5_energy.pdf 
 

# Folder Organization
* `/CRD_File_py_src/` protein internal energy calculation written in python using the professor provided CRD files containing pre-fetched atom variables
* `/PDB_File_py_src/` protein internal energy calculation written in python using self-provided atom variables
* `/data/` contains pdb and crd files 
* `/Share/` contains the Written Report (PDF and DOCX) and Rough Interface along with any supplemental files.
* `/cpp_src/` protein internal energy calculation written in C++


#PLEASE NOTE:- Depending on where you save the files you will have to change the file path in main.py in the PDB_File_py_src folder


# Run time tests

Run with python:

# `/PDB_File_py_src/` Environment Set Up
## Windows
download https://bootstrap.pypa.io/get-pip.py
```shell
python get-pip.py
```

## MSYS2
```shell
pacman -S python3-pip
```

## MacOS
```shell
sudo easy_install pip
sudo pip install --upgrade pip
```

## Get Libraries
```shell
pip install biopython
```

# `/PDB_File_py_src/` Prints each atom in the protein
Use python 3
```shell
cd PDB_File_py_src
python main.py
```


# Run Program to Compute Structure Energy/Metric/Score
Use python 3
```shell
cd CRD_File_py_src
python main.py
```

Please use main.py in CRD_File_py_src to obtain the energy for both conformations of the protein


#!/bin/bash

conda install python==3.7
conda install ipykernel
conda install xarray zarr sparse
conda install matplotlib
conda install plotnine==0.8.0
conda install dask==2021.10.0
conda install h5py
conda install scanpy==1.8.2

pip install --ignore-installe cffi
pip install rpy2

conda install radian
conda install psycopg2
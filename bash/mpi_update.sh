#!/bin/bash

cd  logs/mpi_exp
bash logs.sh
cd ../..
cd  jobs/mpi_exp
bash jobs.sh
cd ../..
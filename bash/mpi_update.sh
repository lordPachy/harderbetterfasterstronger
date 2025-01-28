#!/bin/bash

cd  logs/mpi_bcast
bash logs.sh
cd ../..
cd  jobs/mpi_bcast
bash jobs.sh
cd ../..
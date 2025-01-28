#!/bin/bash

mv pthreads/jobs.sh jobs/pthreads
mv pthreads/logs.sh logs/pthreads
mv sequential/jobs.sh jobs/sequential
mv sequential/logs.sh logs/sequential
mv mpi/jobs.sh jobs/mpi
mv mpi/logs.sh logs/mpi
mv mpi_exp/jobs.sh jobs/mpi_exp
mv mpi_exp/logs.sh logs/mpi_exp
mv mpi_bcast/jobs.sh jobs/mpi_bcast
mv mpi_bcast/logs.sh logs/mpi_bcast
mv pthreads/runs* ./
mv sequential/runs* ./
mv mpi/runs* ./
mv mpi_bcast/runs* ./
mv mpi_exp/runs* ./



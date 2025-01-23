#!/bin/bash

mv pthreads/jobs.sh jobs/pthreads
mv pthreads/logs.sh logs/pthreads
mv sequential/jobs.sh jobs/sequential
mv sequential/logs.sh logs/sequential
mv mpi/jobs.sh jobs/mpi
mv mpi/logs.sh logs/mpi
mv pthreads/runs* ./
mv sequential/runs* ./
mv mpi/runs* ./



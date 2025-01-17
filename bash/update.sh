#!/bin/bash

cd  logs/pthreads
bash logs.sh
cd ..
cd  sequential
bash logs.sh
cd ../..
cd  jobs/pthreads
bash jobs.sh
cd ..
cd  sequential
bash jobs.sh
cd ../..
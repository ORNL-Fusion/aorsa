#!/bin/bash

mpirun -n 1 ../../xaorsa2d &> logfile
diff out38 gold-out38
diff out138 gold-out138
diff rho gold-rho


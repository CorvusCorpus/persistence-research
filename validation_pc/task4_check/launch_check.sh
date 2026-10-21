#!/bin/bash
# Detached launcher for the CPU second pass (run_shells.sh cpu) so it outlives the 2 h background-task limit.
# usage: launch_check.sh "N0:N ..." TAG
cd /home/user/persistence-research/work/gsieve || exit 1
setsid nohup ./run_shells.sh /root/scratch/shells "$1" cpu > /root/scratch/check_$2.out 2>&1 < /dev/null &
echo $! > /root/scratch/check_$2.pid
echo "launched pid $(cat /root/scratch/check_$2.pid)"

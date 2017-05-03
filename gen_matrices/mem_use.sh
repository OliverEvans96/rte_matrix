# File Name: mem_use.sh
# Description: Track memory usage
# Author: Oliver Evans
# Created: Mon May 01, 2017 | 10:52pm EDT
# Last Modified: Mon May 01, 2017 | 11:19pm EDT

# Column headers
printf '%10s %10s %10s\n' 'time' 'memory' 'swap'
# Infinite loop - cancel w/ ctrl-C
while true;
do
	# Print available memory (kB)
	free               | 
	# Ignore headers from `free`
	tail -n +2         | 
	# Take 4th column
	awk '{print $4}'   |
 	# Print as column & prepend unix time
	# (seconds since 1/1/1970)	
	xargs printf '%10d %10d %10d\n' $(date +%s) |\
	# Write to stdout & mem.txt simultaneously
	tee mem.txt;
	# Short delay
	sleep 1
done

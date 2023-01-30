#!/bin/bash

cores=20

MAX_X=460
MAX_Y=310

for ((x=0; x<$MAX_X; x++)); do
  for ((y=0; y<$MAX_Y; y++)); do
    ./Compute_SF $x $y 
	  while :; do
	    background=( $(jobs -p))
	    if (( ${#background[@]} < cores )); then
		    break
	    fi
	    sleep 1
	  done
  done
done


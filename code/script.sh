#!/bin/bash

for NM in 500 1000 2000 4000 8000 16000 32000;do
  # ./simple_matrix.x $NM | grep "time" | cut -d : -f 2 | sed 's/\ //g' >> simple_matrix.txt
  ./operator.x $NM | grep "time" | cut -d : -f 2 | sed 's/\ //g' >> operator.txt
done

#EOF

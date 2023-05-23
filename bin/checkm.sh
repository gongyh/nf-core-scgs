#!/bin/bash

if [ $# -ne 6 ]; then
  echo "Usage: call_checkm.sh ./spades fasta lineage_wf checkm_out 16 checkm_txt"
  exit 0
fi

if [ "$3" != "lineage_wf" ]; then
  checkm taxonomy_wf -t $5 -f $6 -x $2 genus $3 $1 $4
else
  checkm lineage_wf -t $5 -r -f $6 -x $2 $1 $4
fi

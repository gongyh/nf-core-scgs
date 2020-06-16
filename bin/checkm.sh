#!/bin/bash

if [ $# -ne 7 ]; then
  echo "Usage: call_checkm.sh ./spades fasta lineage_wf checkm_out 16 checkm_database checkm_txt"
  exit 0
fi

checkm data setRoot $6
if [ "$3" != "lineage_wf" ]; then
  checkm taxonomy_wf -t $5 -f $7 -x $2 genus $3 $1 $4
else
  checkm lineage_wf -t $5 -r -f $7 -x $2 $1 $4
fi


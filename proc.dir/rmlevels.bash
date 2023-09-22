#!/bin/bash
# removes all LEVEL lines from a CARDS file (including iron level comment line)
# Andreas Sander, last update 15 Oct 2021

if [ -e CARDS ] ; then
  cp CARDS CARDS_LVLBAK
  sed '/LEVEL/d; /* keep old iron/d; /* The following/d; /* and can ONLY/d; /* have been used FEDAT/d' CARDS > CARDS_TMP
  mv CARDS_TMP CARDS
else
  echo "Error: CARDS file could not be found"
fi

#!/bin/bash
source ~/.powrconfig || exit 

kn=$1

printf "Are you sure that you want to remove chain${kn}? y/[n]\n"
read
[ "$REPLY" = y ] || exit

printf "Deleting ...\n"

printf "wrdata${kn}\n"
rm -rf $POWR_WORK/wrdata${kn}

printf "scratch/*/${kn}\n"
rm -rf $POWR_WORK/scratch/formal${kn}
rm -rf $POWR_WORK/scratch/modify${kn}
rm -rf $POWR_WORK/scratch/njn${kn}
rm -rf $POWR_WORK/scratch/steal${kn}
rm -rf $POWR_WORK/scratch/wrstart${kn}
rm -rf $POWR_WORK/scratch/wruniq${kn}

printf "output/*/${kn}\n"
rm -rf $POWR_WORK/output/formal${kn}.cpr
rm -rf $POWR_WORK/output/formal${kn}.log
rm -rf $POWR_WORK/output/formal${kn}.out
rm -rf $POWR_WORK/output/formal${kn}.plot

rm -rf $POWR_WORK/output/njn${kn}.cpr
rm -rf $POWR_WORK/output/njn${kn}.out

rm -rf $POWR_WORK/output/set_repeat${kn}.cpr
rm -rf $POWR_WORK/output/set_repeat${kn}.log

rm -rf $POWR_WORK/output/steal${kn}_backup.cpr
rm -rf $POWR_WORK/output/steal${kn}_backup.log

rm -rf $POWR_WORK/output/steal${kn}.cpr
rm -rf $POWR_WORK/output/steal${kn}.log
rm -rf $POWR_WORK/output/steal${kn}.out
rm -rf $POWR_WORK/output/steal${kn}.plot

rm -rf $POWR_WORK/output/wrstart${kn}.cpr
rm -rf $POWR_WORK/output/wrstart${kn}.log
rm -rf $POWR_WORK/output/wrstart${kn}.out
rm -rf $POWR_WORK/output/wrstart${kn}.plot

rm -rf $POWR_WORK/output/wruniq${kn}.cpr
rm -rf $POWR_WORK/output/wruniq${kn}.log
rm -rf $POWR_WORK/output/wruniq${kn}.out
rm -rf $POWR_WORK/output/wruniq${kn}.plot

printf "wrjobs/*/${kn}\n"

rm -f $POWR_WORK/wrjobs/newformal_cards${kn}
rm -f $POWR_WORK/wrjobs/newdatom${kn}
rm -f $POWR_WORK/wrjobs/formal${kn}
rm -f $POWR_WORK/wrjobs/modify${kn}
rm -f $POWR_WORK/wrjobs/njn${kn}
rm -f $POWR_WORK/wrjobs/set_repeat${kn}
rm -f $POWR_WORK/wrjobs/steal${kn}
rm -f $POWR_WORK/wrjobs/steal${kn}_backup
rm -f $POWR_WORK/wrjobs/wrstart${kn}
rm -f $POWR_WORK/wrjobs/wruniq${kn}


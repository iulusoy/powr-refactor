#!/bin/bash

source ~/.powrconfig || exit

# copy Model data from recent directory to the e.g. wrdata1 directory
# usage: loadmod2.com 1

NR=$1
cp -fp CARDS        $POWR_WORK/wrdata$NR/CARDS 
cp -p xFORMAL_CARDS  $POWR_WORK/wrdata$NR/FORMAL_CARDS 
cp -fp MODEL        $POWR_WORK/wrdata$NR/MODEL 
cp -p DATOM         $POWR_WORK/wrdata$NR/DATOM 
cp -fp FGRID        $POWR_WORK/wrdata$NR/FGRID
cp -pR FEDAT        $POWR_WORK/wrdata$NR/FEDAT
cp -pR FEDAT_FORMAL $POWR_WORK/wrdata$NR/FEDAT_FORMAL
# 
chmod u+w           $POWR_WORK/wrdata$NR/CARDS
chmod u+w           $POWR_WORK/wrdata$NR/MODEL
chmod u+w           $POWR_WORK/wrdata$NR/DATOM
chmod u+w           $POWR_WORK/wrdata$NR/FORMAL_CARDS
chmod u+w           $POWR_WORK/wrdata$NR/FGRID
WD=${PWD##*/}
touch  $POWR_WORK/wrdata${NR}/archive.${WD}
echo "Loaded from ${PWD} at" `date` >>  $POWR_WORK/wrdata${NR}/wrdata${tchain}.hist
echo "Files copied to  wrdata$NR and renamed to archive.file"

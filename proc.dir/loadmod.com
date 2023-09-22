#!/bin/bash

source ~/.powrconfig || exit

# copy Model data from recent directory to the e.g. wrdata1 directory
# usage: loadmod.com 1
NR=$1
cp -fp CARDS        $POWR_WORK/wrdata$NR/archive.CARDS 
cp -p FORMAL_CARDS  $POWR_WORK/wrdata$NR/archive.FORMAL_CARDS 
cp -fp MODEL        $POWR_WORK/wrdata$NR/archive.MODEL 
cp -p DATOM         $POWR_WORK/wrdata$NR/archive.DATOM 
cp -fp FGRID        $POWR_WORK/wrdata$NR/archive.FGRID
cp -pR FEDAT        $POWR_WORK/wrdata$NR/archive.FEDAT
cp -pR FEDAT_FORMAL $POWR_WORK/wrdata$NR/archive.FEDAT_FORMAL
cp -fp NEWDATOM_INPUT $POWR_WORK/wrdata$NR/archive.NEWDATOM_INPUT
cp -fp NEWFORMAL_CARDS_INPUT $POWR_WORK/wrdata$NR/archive.NEWFORMAL_CARDS_INPUT
# 
chmod u+w           $POWR_WORK/wrdata$NR/archive.CARDS
chmod u+w           $POWR_WORK/wrdata$NR/archive.MODEL
chmod u+w           $POWR_WORK/wrdata$NR/archive.DATOM
chmod u+w           $POWR_WORK/wrdata$NR/archive.FORMAL_CARDS
chmod u+w           $POWR_WORK/wrdata$NR/archive.FGRID
chmod u+w           $POWR_WORK/wrdata$NR/archive.NEWDATOM_INPUT
chmod u+w           $POWR_WORK/wrdata$NR/archive.NEWFORMAL_CARDS_INPUT
WD=${PWD##*/}
touch  $POWR_WORK/wrdata${NR}/archive.${WD}
echo "Loaded from ${PWD} at" `date` >>  $POWR_WORK/wrdata${NR}/wrdata${tchain}.hist
echo "Files copied to  wrdata$NR and renamed to archive.file"

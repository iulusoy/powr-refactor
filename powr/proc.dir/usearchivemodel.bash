#!/bin/bash
#copy archive model files to main files 
# (not yet doing links automatically)
cp -f archive.CARDS CARDS
cp -f archive.FORMAL_CARDS FORMAL_CARDS
cp -f archive.DATOM DATOM
cp -f archive.MODEL MODEL
cp -f archive.FGRID FGRID
diff archive.FEDAT FEDAT
diff archive.FEDAT_FORMAL FEDAT_FORMAL
echo "Copied archive.* files as standard files."

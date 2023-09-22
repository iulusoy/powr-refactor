#!/bin/bash
# adds all lines from current LEVELCARDS file to a declared
# Levelcards section of a CARDS file. The LEVELCARDS file is
# created with wrh's levelcards script, which can be called
# via 'levelcards' Andreas Sander, last update 15 Oct 2021

if [ -e CARDS ] ; then
  cp CARDS CARDS_OLDLVL
  # Check for levelcards section (case insensitive)
  seccheck=$(grep -i '^---.*Level.*---' CARDS)
  if [ ${#seccheck} -gt 0 ] ; then
    # Case 1: Levelcards section exists
    sed '/^---.*Levelcards.*---/Ia LVLCMARK' CARDS | sed -e '/LVLCMARK/r LEVELCARDS' -e '/LVLCMARK/d' > CARDS_TMP
  else
    # Case 2: No Levelcards section found: add section at the end of file
    cp CARDS CARDS_TMP
    echo '' >> CARDS_TMP 
    echo '--------------------- Levelcards ---------------------------------------' >> CARDS_TMP
    cat LEVELCARDS >> CARDS_TMP
    echo '' >> CARDS_TMP 
  fi
  mv CARDS_TMP CARDS
else
  echo "Error: CARDS file could not be found"
fi

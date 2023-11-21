      PROGRAM MAINnewdatom 
C***  Provide Link data for possible use in the programm
      CHARACTER LINK_DATE*30, LINK_USER*10, LINK_HOST*60
      COMMON / COM_LINKINFO / LINK_DATE, LINK_USER, LINK_HOST
      LINK_DATE = 'Di 21. Nov 13:00:28 CET 2023'
      LINK_USER = 'inga'
      LINK_HOST = 'ssc-laptop01'
                               
      CALL newdatom 
      END
      SUBROUTINE FINDELEMENT (SYMBOL, NAME, STAGE, ATOMICMASS)
C*****************************************************************
C***  Hier wird einem Elementsymbol der Elementname, Masse und Stage
C***  zugeordnet.
C*****************************************************************

      CHARACTER SYMBOL*(*)
      PARAMETER (MAXELEM = 26)
      CHARACTER*10 ELEMNAME(MAXELEM)
      CHARACTER*2 ELEMSYMBOL(MAXELEM)
      REAL ELEMMASSE(MAXELEM)
      REAL ELEMSTAGE(MAXELEM)
C      CHARACTER*5 ELEMMASSE(MAXELEM)
C      CHARACTER*3 ELEMSTAGE(MAXELEM)
      INTEGER Z
      CHARACTER*10 NAME

      
      DATA ELEMSYMBOL /'H ', 'HE', 'LI', 
     2                 'BE', 'B ', 'C ', 
     3                 'N ', 'O ', 'F ',
     4                 'NE', 'NA', 'MG', 
     5                 'AL', 'SI', 'P ',
     6                 'S ', 'CL', 'AR', 
     7                 'K ', 'CA', 'SC', 
     8                 'TI', 'V ', 'CR', 
     9                 'MN', 'G '/

      
      DATA ELEMNAME /'HYDROGEN  ', 'HELIUM    ', 'LITHIUM   ', 
     2               'BERYLLIUM ', 'BORON     ', 'CARBON    ', 
     3               'NITROGEN  ', 'OXYGEN    ', 'FLUORINE  ',
     4               'NEON      ', 'SODIUM    ', 'MAGNESIUM ', 
     5               'ALUMINIUM ', 'SILICON   ', 'PHOSPHORUS',
     6               'SULFUR    ', 'CHLORINE  ', 'ARGON     ', 
     7               'POTASSIUM ', 'CALCIUM   ', 'SCANDIUM  ', 
     8               'TITANIUM  ', 'VANADIUM  ', 'CHROMIUM  ', 
     9               'MANGANESE ', 'GENERIC   '/

C    Fuer H, HE, C, N, O, SI, P, NE sind hier Werte eingetragen, die ich aus
C    DATOM-Files zusammengesucht habe, der Rest ist aus dem
C    Periodensystem abgeschrieben.   
      DATA ELEMMASSE / 1.00, 4.00, 6.94, 
     2               9.01, 10.81, 12.00, 
     3               14.00, 16.00, 19.00,
     4               20.18, 22.99, 24.31, 
     5               26.98, 28.09, 30.97,
     6               32.07, 35.45, 39.95, 
     7               39.10, 40.08, 44.96, 
     8               47.87, 50.94, 52.00, 
     9               54.94, 0.00 /
     
     
C    Fuer H, HE, C, N, O, SI, P, NE sind hier Werte eingetragen, die ich aus
C    DATOM-Files zusammengesucht habe, der Rest ist willkuerlich auf 2 gesetzt     
      DATA ELEMSTAGE / 2. ,3. ,2. , 
     2                 2. ,2. ,5. , 
     3                 6. ,7. ,2. ,
     4                 2. ,2. ,2. , 
     5                 2. ,5. ,5. , 
     6                 2. ,2. ,2. ,  
     7                 2. ,2. ,2. , 
     8                 2. ,2. ,2. , 
     9                 2. ,2. /   
     

      
      DO Z=1, MAXELEM
         IF (SYMBOL .EQ. ELEMSYMBOL(Z)) THEN
           WRITE(NAME,'(A)') ELEMNAME(Z)
C	   READ (ELEMSTAGE(Z), FMT='(F5.0)', ERR=900) STAGE
C	   READ (ELEMMASSE(Z), FMT='(F6.2)', ERR=900) MASSE
	   STAGE = ELEMSTAGE(Z)
	   ATOMICMASS = ELEMMASSE(Z)
           EXIT
         ENDIF
      END DO
      
      RETURN 
      
  900  PRINT *, 'ERROR: ' 
       STOP '*** ERROR DETECTED' 

      END

       FUNCTION IDX(TEXT)
C***  Returns the number of non-blank characters of string TEXT
       CHARACTER*(*) TEXT
       DO IDX=LEN(TEXT),1,-1
          IF(TEXT(IDX:IDX).NE.' ') RETURN
       END DO
       IDX=0
       RETURN
       END

      SUBROUTINE JSYMSET (SYMBOL,WORD)
C***********************************************************************
C***  USER -DEFINED UNICOS VERSION OF THE COS LIBRARY ROUTINE JSYMSET
C***  OPENS A FILE NAMED AS THE SPECIFIED COS SYMBOL AND
C***  WRITES THE CHARACTER STRING "WORD" INTO THIS FILE
C***********************************************************************
      CHARACTER SYMBOL*(*), WORD*(*)

      OPEN (66,FILE=SYMBOL, STATUS='UNKNOWN')
      WRITE (66,'(A)') WORD
      CLOSE (66)

      RETURN
      END
      SUBROUTINE NEWDATOM
C**********************************************************************
C***
C***  Program to create a DATOM (= atomic data) file for the PoWR code
C***      The program needs a database of separate DATOM files 
C***      for each ion to be included. 
C***
C***  This program requires an input file NEWDATOM_INPUT
C***  to specify which elements and ions and how many levels 
C***  are to be included.
C***
C**********************************************************************
      PARAMETER (NDIM=1000) !Maximale Zahl der Levels
      PARAMETER (MAXION=27)
      CHARACTER LEVEL(NDIM, MAXION)*10
      CHARACTER ZEILE*80, NAMEION*5(MAXION), FILENAME*255 
      CHARACTER ELEMENT*2 !Elementsymbol
      CHARACTER IONDEGREE(MAXION)*10 !vorkommende Ionen
      CHARACTER HIGHESTION*10, ROMAN(20)*5
      CHARACTER ACTPAR*20, ACTPARION*20, ACTPAR2*20
      CHARACTER ATOM*80
      CHARACTER DEGREE*20
      LOGICAL DRTRANSITGLOBAL, KSHELLGLOBAL, DRTRANSITELEM, KSHELLELEM
      LOGICAL, DIMENSION(MAXION) :: DRTRANSIT, KSHELL
      CHARACTER*10 NAME !Elementname
      REAL STAGE
      DIMENSION NLEVEL(MAXION)
      CHARACTER ZEILE1*80
      CHARACTER*10 CDATE, CTIME
      REAL MASSE
      DIMENSION LEVELCOUNT(MAXION) 
      CHARACTER PATHSAVE*255
      CHARACTER PATH(MAXION)*255
      LOGICAL BTRUEAUTOLEV
      CHARACTER FEHLERZEILE(MAXION)*80

C***  Link data to identify program version
      CHARACTER LINK_DATE*30, LINK_USER*10, LINK_HOST*60
      COMMON / COM_LINKINFO / LINK_DATE, LINK_USER, LINK_HOST

C***  Converting roman numerals to arabic numerals
      DATA ROMAN / 'I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII',
     >             'IX', 'X', 'XI', 'XII', 'XIII', 'XIV', 'XV', 
     >             'XVI', 'XVII', 'XVIII', 'XIX', 'XX' /  

C***  Write Link Data (Program Version) tp CPR file
      WRITE(0,'(2A)') '>>> NEWDATOM: Program Version from ', 
     >                LINK_DATE
      WRITE(0,'(4A)') '>>> created by ', LINK_USER(:IDX(LINK_USER)),
     >      ' at host ', LINK_HOST(:IDX(LINK_HOST))
    

C***  Error default for PATH to the data
      PATHSAVE = '..UNDEFINED..'


C     Ausgabe in einer neuen Datei 
      OPEN (99,FILE='DATOM', ERR = 904)  
    
C     Einlesen des Steuerfiles
      OPEN (98,FILE = 'NEWDATOM_INPUT', ACTION = 'READ', ERR = 905) 
      
C     Gesamtes Steuerfile als Kommentar oben drueber schreiben, d.h.
C     einlesen und ausgeben
      CALL DATE(CDATE)
      CALL TIME(CTIME)
      WRITE (99,'(A)') 
     > '* This DATOM-FILE has been created at ' // CDATE 
     >                  // '' // CTIME
      WRITE (99,'(A)') '* with the program NEWDATOM' 
      WRITE (99,'(A)') '*'
      WRITE (99,'(A)') '* NEWDATOM_INPUT was'       
    7 DO
       READ (98,'(A)', END=1) ZEILE1
       IF (ZEILE1(:1) .EQ. '*' .OR. ZEILE1(:1) .EQ. '-') GOTO 7
       WRITE (99,'(A)') '* ' // ZEILE1(:IDX(ZEILE1))
      END DO
    1 CONTINUE
      WRITE (99,'(A)') '*'      

      REWIND(98)   ! Am Schluss den Pointer wieder an den Anfang setzen
      
C     Defaults 
      DRTRANSITGLOBAL = .FALSE.
      BTRUEAUTOLEV = .TRUE.
      KSHELLGLOBAL = .TRUE.
      
C     Schleife fuer die Elemente    
      DO 
    3  READ (98,'(A)', END = 5) ZEILE !Datei zeilenweise durchlesen
       IF (ZEILE .EQ. '') GOTO 3 !Leerzeilen abfangen 
       CALL SARGV(ZEILE,1,ACTPAR)   !Erstes Wort der Zeile (Keyword)
C      einlesen und anschliessend untersuchen (Fallunterscheidung)

C      Path-Variable einlesen in PATHSAVE
       IF (ACTPAR .EQ. 'PATH') THEN
         CALL SARGV(ZEILE,2,PATHSAVE) 

C      Hier wird DRTRANSIT global fuer alle Ionen der folgenden 
C       Elemente angeschaltet  
       ELSEIF (ACTPAR .EQ. 'DRTRANSIT') THEN 
        DRTRANSITGLOBAL = .TRUE. 
C       Test for second parameter to include levels below threshold
        CALL SARGC (ZEILE, NARG)
        IF (NARG .GT. 1) THEN
           CALL SARGV(ZEILE,2,ACTPAR2) 
           IF (ACTPAR2(:11) .EQ. 'INCLUDE_NEG') BTRUEAUTOLEV=.FALSE.
        ENDIF

C      Hier wird DRTRANSIT global fuer alle Ionen der folgenden 
C       Elemente abgeschaltet  
       ELSEIF (ACTPAR .EQ. 'NO-DRTRANSIT') THEN 
        DRTRANSITGLOBAL = .FALSE. 

       ELSE IF (ACTPAR .EQ. 'NO-K-SHELL') THEN
        KSHELLGLOBAL = .FALSE.       

       ELSE IF (ACTPAR .EQ. 'K-SHELL') THEN
        KSHELLGLOBAL = .TRUE.       
        
C      Einzelnes Element wird hier drin abgearbeitet
       ELSE IF (ACTPAR .EQ. 'ELEMENT') THEN
          NION = 0 !Zaehlt, wie viele Ionen angegeben sind
          CALL SARGC (ZEILE, NPAR) !Anzahl der Parameter ermitteln
          IF (NPAR .LT. 2) GOTO 900 !Fehler abfangen
          CALL SARGV(ZEILE,2,ELEMENT) !Element = zweites Wort der Zeile
  
C         Ausnahme fuer Generic: Hier wird nur Ionlow und Iontop aus dem
C         Steuerfile eingelesen und hingeschrieben, alles andere ist fest
  
          IF (ELEMENT .EQ. 'G') THEN
           CALL SARGC (ZEILE, NPAR) !Anzahl der Parameter ermitteln
           IF (NPAR .LT. 4) GOTO 903 !Fehler abfangen
           CALL SARGV(ZEILE,3,ACTPAR) !Ionlow steht an dritter Stelle in
C                                 der entsprechenden Zeile im Steuerfile
           READ (ACTPAR,'(I20)') LOW  !Typumwandlung
           CALL SARGV(ZEILE,4,ACTPAR) !Iontop steht an vierter Stelle
           READ (ACTPAR,'(I20)') NUP  !Typumwandlung
C          Ausgabe der Daten fuer Generic
           WRITE (99,'(A)')
     >        '*KEYWORD--  ---NAME--- SYMB   IONLOW  IONTOP'
           WRITE (99,'(A)')
     >        '**************************************************'
           WRITE (99,'(A27,7X,I2,6X,I2)')
     >       'ELEMENT     GENERIC    (G )', LOW, NUP
           WRITE (99,'(A)')
     >       '*           =======                              *'
           WRITE (99,'(A)')
     >       '**************************************************'
   
C         Andere Elemente
          ELSE

C          We have to find the highest ionization stage
           NBACKSPACE = 0
           DO
             READ (98,'(A)',END=41) ZEILE
             IF (ZEILE .EQ. '') GOTO 42
             CALL SARGV (ZEILE,1,ACTPARION)
             IF (ACTPARION .EQ. 'ION') THEN
                CALL SARGC (ZEILE,NPAR) 
                IF (NPAR .LT. 2) GOTO 902 ! Fehler
                CALL SARGV(ZEILE,2,HIGHESTION) 
                NBACKSPACE = NBACKSPACE + 1
             ELSEIF (ACTPARION .EQ. 'ELEMENT') THEN
 41             NBACKSPACE = NBACKSPACE + 1
                EXIT
             ELSE 
 42             NBACKSPACE = NBACKSPACE + 1
             ENDIF             
           ENDDO
C          Move REC-Pointer back to first line of ELEMENT
           DO IBACKSPACE = 1, NBACKSPACE
              BACKSPACE 98
           ENDDO
           DO I=1, 20
              IF (HIGHESTION .EQ. ROMAN(I)) THEN
                 MAINSTAGE = MAX((I-1),1)
                 EXIT
              ENDIF
           ENDDO
              
C          Aufruf einer Subroutine, die aus dem Elementsymbol auf
C          weitere Eigenschaften schliesst  
           CALL FINDELEMENT(ELEMENT, NAME, STAGE, ATOMICMASS)
C          Ausgabe des Headers fuer das aktuelle Element
           WRITE (99,'(A)')
     >       '*======================================================='
           WRITE (99,'(A)')
     >       '*KEYWORD--  ---NAME--- SYMB   ATMASS   STAGE' 
           WRITE (99,'(A)')
     >       '************************************************'
C              USE MAINSTAGE INSTEAD OF STAGE
           WRITE (99, '(A7,5X,A10,1X,A1,A2,A1,3X,F6.2,3X,I3)') 
     >       'ELEMENT', NAME, '(', ELEMENT, ')', ATOMICMASS, MAINSTAGE
  
           WRITE (99,'(A)')
     >       '*           ======                             *'
           WRITE (99,'(A)')
     >       '************************************************'

C          default: global settings for DR-TRANSIT and K-SHELL in element
           DRTRANSITELEM = DRTRANSITGLOBAL
           KSHELLELEM = KSHELLGLOBAL
           
C***       Read ION options for each ion
            
    4      READ (98,'(A)', END = 2) ZEILE !Ionenzeilen nacheinander
C                                          einlesen
           IF (ZEILE .EQ. '') GOTO 4 !Leerzeilen abfangen
           CALL SARGV(ZEILE,1,ACTPARION) !erstes Wort ist Schluesselwort
   
C          Pfadvariable aendert sich
           IF (ACTPARION .EQ. 'PATH') THEN
C                              ====
            CALL SARGV(ZEILE,2,PATHSAVE) !PATH neu setzen
    
C          Kommandos zwischen ELEMENT und vor einem ION-Kommando 
C           koennen DR-TRANSIT und K-SHELL fuer alle folgenden Ionen
C           desselben Elements ein- oder ausschalten
           ELSEIF (ACTPARION == 'DRTRANSIT') THEN 
             DRTRANSITELEM = .TRUE. 

           ELSEIF (ACTPARION == 'NO-DRTRANSIT') THEN 
             DRTRANSITELEM = .FALSE. 

           ELSEIF (ACTPARION == 'NO-K-SHELL') THEN
             KSHELLELEM = .FALSE.       

           ELSEIF (ACTPARION == 'K-SHELL') THEN
             KSHELLELEM = .TRUE.       
    
C          Ionenzeile, hier werden Ionname und Parameter fuer einzelne 
C          Ionen angegeben
           ELSE IF (ACTPARION .EQ. 'ION') THEN
C                                   ===
             CALL SARGC (ZEILE, NPAR) !Zaehle die Parameter der Zeile
             IF (NPAR .LT. 2) GOTO 902 !Fehler abfangen
             NION = NION + 1 !Ionenzaehler hochsetzen fuer Schleifen in
C                             newdatomion
             FEHLERZEILE(NION) = ZEILE

             IF (PATHSAVE .EQ. '..UNDEFINED..') GOTO 908

             PATH(NION) = PATHSAVE !aktueller Pfad           

             CALL SARGV(ZEILE,2,IONDEGREE(NION)) !IONDEGREE einlesen     

C            entscheiden, ob DRTRANSIT mitgenommen wird (default = element setting)
             DRTRANSIT(NION) = DRTRANSITELEM

C            entscheiden, ob K-SHELL mitgenommen wird (default = element setting)
             KSHELL(NION) = KSHELLELEM

C            Zahl der Levels einschraenken oder (wenn keine
C            Einschraenkung gegeben ist) alle nehmen
             NLEVEL(NION) = -1 !Default: Wenn kein Level angegeben
C                               ist, nimmt er alle     
             DO ICOUNT = 1, NPAR !Gehe alle Parameter durch und pruefe,
C                                 ob irgendwo DRTRANSIT, NO-DRTRANSIT
C                                 oder NLEVEL steht      
              CALL SARGV(ZEILE,ICOUNT,ACTPAR) !Zeile wird in Parameter
C                                              zerlegt
              IF (ACTPAR .EQ. 'NO-DRTRANSIT') THEN
               DRTRANSIT(NION) = .FALSE. !DRTRANSIT fuer das ION
C                                         ausschalten
              ELSE IF (ACTPAR .EQ. 'DRTRANSIT') THEN
               DRTRANSIT(NION) = .TRUE. !DRTRANSIT fuer dieses ION
C                                        einschalten       
              ELSE IF (ACTPAR .EQ. 'NO-K-SHELL') THEN
               KSHELL(NION) = .FALSE.   !K-SHELL fuer das ION
C                                         ausschalten
              ELSE IF (ACTPAR .EQ. 'K-SHELL') THEN
               KSHELL(NION) = .TRUE.    !K-SHELL fuer dieses ION
C                                        einschalten       
              ELSE IF (ACTPAR .EQ. 'NLEVEL') THEN !Zahl der Level
C                                                  beschraenken
               IF (ICOUNT+1 .GT. NPAR) GOTO 906
               CALL SARGV(ZEILE,ICOUNT+1,ACTPAR) !Levelzahl wird
C                           fuer jedes Ion als Character abgespeichert
               READ (ACTPAR, '(I20)', ERR = 907) NLEVEL(NION)

              END IF 
             END DO      
           
           ELSEIF (ACTPARION .EQ. 'ELEMENT') THEN 
C             Reading for one Element is finished: 
C             Reset the pointer for continuing later with next element
              BACKSPACE (98)     
              GOTO 2
     
           ENDIF  

C***       Read next line within ELEMENT block
           GOTO 4 

C      Reading input options for last ELEMENT completed
C      (next ELEMENT or E-O-F encountered) 
C      Subroutine NEWDATOMION wird fuer jedes Element ein Mal
C      aufgerufen, Informationen fuer einzelne Ionen stehen in Arrays 
C      und werden hier uebergeben   

    2      CALL NEWDATOMION (ELEMENT,NION,IONDEGREE,
     >          PATH,DRTRANSIT,NLEVEL,LEVEL,NDIM,MAXION, 
     >          LEVELCOUNT, FEHLERZEILE, KSHELL, BTRUEAUTOLEV)
         
          END IF !Element finished
                 
       ENDIF  
               
      END DO !End of loop reading input lines
    5 CONTINUE 
      
      CLOSE (99) !Schliessen der output-DATOM-Datei
      CLOSE (98) !Schliessen des Steuerfiles

C**   Regular program end - remove error code 
   99 CONTINUE
      CALL JSYMSET ('G0', '0')
      STOP 'O.K.'

C********  ERROR EXITS **************************************
  900 WRITE(0,*) 'ERROR: ELEMENT option without parameter'
      GOTO 990
      
  901 WRITE(0,*) 'ERROR: ELEMENT GENERIC NOT FOUND'
      GOTO 990    
      
  902 WRITE(0,*) 'ERROR: ION option without parameter'
      GOTO 990     
  
  903 WRITE(0,*) 'ERROR: ELEMENT GENERIC without parameter'
      GOTO 990         
      
  904 WRITE(0,*) 'ERROR: CANNOT WRITE IN FILE DATOM'
      GOTO 990     
      
  905 WRITE(0,*) 'ERROR: NEWDATOM_INPUT NOT FOUND'
      GOTO 990       

  906 WRITE(0,*) 'ERROR: LEVEL option without parameter'
      GOTO 990   
  
  907 WRITE(0,*) 'ERROR: invalid parameter for NLEVEL'
      GOTO 990      
     
  908 WRITE(0,*) 'ERROR: PATH to atomic data not (yet) defined'
      GOTO 990
          
            
  990 WRITE(0,*) '*** The error occured in the following line:'
      WRITE(0,*) ZEILE(:IDX(ZEILE))
      STOP '*** ERROR DETECTED IN PROGRAM NEWDATOM'
      
      
      END
      SUBROUTINE NEWDATOMION (ELEMENT,NION,IONDEGREE,PATH,
     >           DRTRANSIT,NLEVEL,LEVEL,NDIM, MAXION, LEVELCOUNT,
     >           FEHLERZEILE, KSHELL, BTRUEAUTOLEV)
C**********************************************************************
C***
C***  Subroutine called by NEWDATOM for each element which is to be included, 
C***  reads data from separate DATOM files for each ion
C***  and writes it into the new DATOM output-file in the following order:
C***  LEVELS, LINES, CONTINUUM, DRTRANSIT, K-SHELL 
C***  
C**********************************************************************     

      CHARACTER(255) :: FILENAME
      CHARACTER(80) :: ZEILE, ZEILE2, ZEILE3, ZEILE4
      CHARACTER(5), DIMENSION(20) :: NAMEION 
      CHARACTER(2) :: ELEMENT
      INTEGER, DIMENSION(MAXION) :: NLEVEL, LEVELCOUNT
      LOGICAL, DIMENSION(MAXION) :: DRTRANSIT, KSHELL
      CHARACTER LEVEL(NDIM, MAXION)*(*)
      LOGICAL :: LEVELFOUND, LEVELFOUNDA, LEVELFOUNDZ, ZIEL
      CHARACTER(255), DIMENSION(MAXION) :: PATH
      CHARACTER(80), DIMENSION(MAXION) :: FEHLERZEILE
      CHARACTER(10), DIMENSION(MAXION) :: IONDEGREE
      CHARACTER(10), DIMENSION(MAXION, NDIM) :: ZIELNIVEAU
      INTEGER :: Z
      LOGICAL :: BTRUEAUTOLEV



C     Oeffnen aller Ionendateien fuer ein Element
      
   1  DO ION=1, NION    
         FILENAME=PATH(ION)(:IDX(PATH(ION))) // '/' // 'DATOM.' // 
     >   ELEMENT(:IDX(ELEMENT)) // '_' // IONDEGREE(ION) 
         OPEN (ION,FILE=FILENAME, ACTION = 'READ', ERR=901)
         WRITE (0,'(2A)') 'Opened: ', FILENAME(:IDX(FILENAME))
      ENDDO

      
C     Um als Zielniveaus benoetigte Level vorher zu finden: 
C     Vorbelegen des Arrays mit Leerstrings
      DO ION=1, NION
       DO I=1,NDIM
        ZIELNIVEAU(ION-1,I) = ''
       END DO
      END DO      
      
C     Um vorher zu wissen, welche Level ich zusaetzlich brauche, muss
C     ich schon vorher die Continuumskarten lesen und alle benoetigten
C     Zielniveaus fuer jedes Ion und jedes Level zwischenspeichern.     

      DO ION=1, NION
       Z=0
       DO
        READ (ION,'(A)', END=8 ) ZEILE2
        IF (ZEILE2(:10) .EQ. 'CONTINUUM ') THEN
         Z=Z+1  !Nummer der aktuellen Levelzeile innerhalb eines Ions  
         ZIELNIVEAU(ION,Z)=ZEILE2(70:79) !Zielniveau abspeichern
        END IF       
       END DO
    8  REWIND(ION) !Pointer wieder zurueck auf den Anfang setzen
      END DO
      
C     Hier ist der eigentliche Anfang      
     
C     Level
      WRITE (99,'(A)') !Header fuer Levelkarten ausgeben
     >      '*KEYWORD--  ---NAME--- CH WEIG--ENERGY-- EION----- QN'
      DO ION=1, NION !Schleife ueber alle Ionen
C        WRITE(0,*) 'Levelmax =', NLEVEL(ION)    !Testausgabe   
        LEVELCOUNT(ION) = 0 !Zaehler fuer die Anzahl der Levels
       DO  !Schleife fuer jedes einzelne Ion       
         READ (ION,'(A)', END=2 ) ZEILE2 !Fuer jedes Ion seine
C                                         Datendatei lesen
         IF (ZEILE2(:10) .EQ. 'LEVEL     ') THEN !Levelzeilen raussuchen
          IF (LEVELCOUNT(ION)+1 .GT. NDIM) GOTO 900 !Fehler durch zu
C                                                   viele Levels abfangen
           
C          Pruefen, ob noch Levels als Ziel benoetigt werden oder ob das
C          Array mit den Zielniveaus schon leer geraeumt wurde
           ZIEL = .TRUE.
           DO I=1,NDIM
            ZIEL = ZIEL .AND. ZIELNIVEAU(ION-1,I) .EQ. ''
           END DO 
   
C         Fuer das erste Ion gibt es keine niedrigeren, deshalb prueft 
C         man nur, ob die gewuenschte Zahl an Levels abgearbeitet wurde:
          IF (ION .EQ. 1) THEN 
           IF (LEVELCOUNT(ION)+1 .GT. NLEVEL(ION) .AND. 
     >      NLEVEL(ION) .NE. -1) GOTO 2

C         Fuer alle anderen muss man pruefen, ob die Levelzahl erreicht
C         ist und alle benoetigen Levels aus dem niedrigeren Ion
C         geschrieben wurden
          ELSE    
            IF (LEVELCOUNT(ION)+1 .GT. NLEVEL(ION) .AND. 
     >       NLEVEL(ION) .NE. -1 .AND. ZIEL) GOTO 2
          END IF
   
C          Mitzaehlen, wie viele Levels geschrieben werden:
           LEVELCOUNT(ION) = LEVELCOUNT(ION) + 1 
C          Merken, welches Level gleich geschrieben wird
C           (wird fuer Lines und Continuum benoetigt):
           LEVEL (LEVELCOUNT(ION),ION) = ZEILE2(13:22) 
           WRITE (99,'(A)') ZEILE2(:IDX(ZEILE2)) !Levelzeile

C          Wenn ein Zielniveau geschrieben wurde, wird es aus dem Array
C          geloescht und durch einen Leerstring ersetzt
           DO J=1, NDIM
            IF (ZEILE2(13:22) .EQ. ZIELNIVEAU (ION-1,J))THEN
             ZIELNIVEAU (ION-1,J) = '' 
            END IF   
           END DO

         ELSE IF (ZEILE2(:10) .EQ. 'LINE      ') THEN 
C          Jetzt sind die Levels zuende und es kommen Lines
           IF (LEVELCOUNT(ION) .LT. NLEVEL(ION)) THEN
             WRITE(0,*) 'WARNING: MORE LEVELS REQUESTED THAN AVAILABLE'
             WRITE(0,'(5A,I3,A,I3)') 
     >            '  Ion: ', ELEMENT,' ', IONDEGREE(ION),  
     >            '  Requested: ', NLEVEL(ION), 
     >            '  Available: ', LEVELCOUNT(ION)
           ENDIF
 
           BACKSPACE (ION) !Pointer einen zurueck um auch die erste
C                           Line-Zeile zu erwischen
           EXIT !Wenn die Linien anfangen, soll er hier die
C                Level-Schleife verlassen
         
         ELSE IF (ZEILE2(:10) .EQ. 'CONTINUUM ') THEN !Falls keine
                              !Linien da sind, kommt jetzt das Continuum
           BACKSPACE (ION) !Pointer einen zurueck um auch die erste
C                           Continuum-Zeile zu erwischen
           EXIT !Wenn es keine Linien mehr gibt und nach den Levels das
C           Continuum anfaengt, soll er hier die Line-Schleife verlassen
  
         ENDIF !Ende der Unterscheidung zwischen Levels, Lines,
C               Continuum
          
         CYCLE !Weitere Zeilen fuer ein Ion einlesen
    2    EXIT !Ende der Datei erreicht         
       END DO !Ende der Schleife fuer jedes einzelne Ion

      END DO !Ende der Schleife ueber alle Ionen
      
C     Line 
      WRITE (99,'(A)') !Header fuer Line-Karten ausgeben
     >'*KEYWORD--UPPERLEVEL  LOWERLEVEL--EINSTEIN  RUD-CEY ''
     >--COLLISIONAL COEFFICIENTS--'
      DO ION=1, NION !Schleife ueber alle Ionen
       DO !Einlesen der Datei fuer jedes einzelne Ion
        READ (ION,'(A)',END=3) ZEILE2 
        IF (ZEILE2(:10) .EQ. 'LINE      ') THEN  !Line-Zeilen

C     Ausgangs- und Zielniveau muessen existieren
         LEVELFOUNDA = .FALSE.
         LEVELFOUNDZ = .FALSE.
      
         DO I=1, LEVELCOUNT(ION) !Nach Ausgangs- und Zielniveau suchen
          LEVELFOUNDA = LEVELFOUNDA .OR. LEVEL(I,ION) .EQ. ZEILE2(23:32)
          LEVELFOUNDZ = LEVELFOUNDZ .OR. LEVEL(I,ION) .EQ. ZEILE2(11:20)  
         END DO
      
         IF (LEVELFOUNDA .AND. LEVELFOUNDZ) THEN !Nur wenn beide
C                                                 gefunden wurden, schreibe die Line-Zeile
         WRITE (99,'(A)') ZEILE2(:IDX(ZEILE2))
         END IF
        ELSE IF (ZEILE2(:10) .EQ. 'CONTINUUM ') THEN !Lines sind vorbei,
C                                                     Continuum faengt an
          BACKSPACE (ION)
          EXIT ! Wenn er eine Continuumszeile erreicht, soll er die
C                Lines-Schleife verlassen
        ENDIF !Ende der Fallunterscheidung zw. Lines und Continuum
        CYCLE !gleiche Datei, naechste Zeile
    3  EXIT !Ende der Datei erreicht
       END DO !Ende der Schleife fuer jedes einzelne Ion
      END DO !Ende der Schleife ueber alle Ionen
      
C     Continuum           
      WRITE (99,'(A)') '*KEYWORD  LOWERLEVEL ----SIGMA ----ALPHA
     >----SEXPO -IGAUNT- -KEYCBF- --IONLEV--' !Header schreiben
      DO ION=1, NION-1 !Schleife ueber Ionen, aber beim letzten Ion 
C      duerfen keine Kontinuumskarten mehr geschrieben werden, 
C      weil dann die Zielniveaus fehlen.
       Z=0
       DO !Schleife ueber die Zeilen jeder einzelnen Iondatei
        READ (ION,'(A)',END=4) ZEILE2 !Datei fuer das Ion einlesen
        IF (ZEILE2(:10) .EQ. 'CONTINUUM ') THEN !Continuumszeilen

C     Folgezeilen mitnehmen   
         IF (ZEILE2(52:57) .EQ. 'PIKB12' .OR. ZEILE2(52:59) .EQ. 'OPAPROIX') THEN
          READ (ION,'(A)',END=4) ZEILE3 !naechste Zeile zwischenspeichern
          BACKSPACE(ION) !Cursor wieder zurueck
         END IF 

C     Ausgangslevel muss existieren und gewollt sein  
        LEVELFOUND = .FALSE.
         DO I=1, LEVELCOUNT(ION) !Sucht nach dem Ausgangslevel
          LEVELFOUND = LEVELFOUND .OR. LEVEL(I,ION) .EQ. ZEILE2(11:20)  
         END DO
 
         IF (LEVELFOUND) THEN 
          WRITE (99,'(A)') ZEILE2(:IDX(ZEILE2)) !Schreibe die Zeile
  
          IF (ZEILE2(52:57) .EQ. 'PIKB12' .OR. 
     >     ZEILE2(52:59) .EQ. 'OPAPROIX') THEN
             WRITE (99,'(A)') ZEILE3(:IDX(ZEILE3)) !Schreibe die Folgezeile
          END IF
     
         END IF 
        ELSE IF (ZEILE2(:10) .EQ. 'DRTRANSIT ') THEN !Continuum vorbei,
C                                                     DRTRANSIT beginnt
          BACKSPACE (ION) !Cursor zuruecksetzen um keine Zeile
                          !auszulassen
          EXIT !Continuumsschleife verlassen
        ELSE IF (ZEILE2(:10) .EQ. 'K-SHELL   ') THEN !Kein DRTRANSIT
C                                                     vorhanden, also gleich zu K-SHELL
          BACKSPACE (ION)
          EXIT !Continuumsschleife verlassen
        ENDIF !Ende der Fallunterscheidung zw. Continuum, DRTRANSIT und
C              K-SHELL
        CYCLE !naechste Zeile
    4   EXIT !Ende der Datei erreicht
       END DO !Ende der Schleife ueber Zeilen fuer einzelnes Ion
      END DO !Ende der Schleife ueber alle Ionen

C     DRTRANSIT
      J=1 !Hilfsvariable, nur fuer die Kopfzeile von Bedeutung
      DO ION=1, NION !Schleife ueber alle Ionen
        IF (DRTRANSIT(ION)) THEN !Nur wenn DRTRANSIT an ist
         DO !Schleife ueber alle Zeilen eines Ions
          READ (ION,'(A)',END=5) ZEILE2 !Einlesen
          IF (ZEILE2(:10) .EQ. 'DRTRANSIT ') THEN !nur DRTRANSIT-Zeilen
C                                                  nehmen  
C     Ausgangslevel muss existieren und gewollt sein
           LEVELFOUND = .FALSE.
           DO I=1, LEVELCOUNT(ION) !Suche nach dem Ausgangslevel
            LEVELFOUND = LEVELFOUND .OR. LEVEL(I,ION) .EQ. ZEILE2(11:20)  
           END DO
           IF (LEVELFOUND) THEN  !Wenn Ausgangslevel gefunden
C           Schreibe beim ersten Mal den Header   
            IF (J .EQ. 1) WRITE (99,'(A)') '*KEYWORD  LOWERLEVEL'
     >       // '  UPPERLEVEL --G- --ENERGY-- -EINSTEIN-'

            J = J+1 !Hochzaehlen, damit der Header nur einmal
C                    geschrieben wird

C***        Copy DRTRANSIT line to output
C***        If BTRUEAUTOLEV=.TRUE. (default), those entries with
C***           negative level energy (relative to ionization threshold)
C***           are skipped (since 25-Jul-2023, wrh)
              READ (ZEILE2, '(38X,F10.0)') EAUTO
              IF (.NOT. BTRUEAUTOLEV .OR. EAUTO .GT. .0) THEN
               WRITE (99,'(A)') ZEILE2(:IDX(ZEILE2)) !Schreibe die Zeile
              ENDIF 
           END IF
    
          ELSE IF (ZEILE2(:10) .EQ. 'K-SHELL   ') THEN !DRTRANSIT
C                                           vorbei, K-SHELL beginnt
            BACKSPACE (ION)
            EXIT !DRTRANSIT-Schleife verlassen
          END IF
          CYCLE !naechste Zeile
    5     EXIT !Ende der Datei erreicht
         END DO !Ende der Schleife fuer einzelnes Ion
       ELSE IF (ZEILE2(:10) .EQ. 'K-SHELL   ') THEN !Auch wenn es kein
C            DRTRANSIT gab, koennte irgendwann ein K-SHELL-Block kommen
           BACKSPACE (ION)
           EXIT !DRTRANSIT-Schleife verlassen
       END IF  !Ende der Fallunterscheidung zwischen DRTRANSIT und
C               K-SHELL 
      END DO !Ende der Schleife ueber alle Ionen
       
C     K-SHELL 
    7 DO ION=1, NION !Fuer jedes Ion seine Datendatei einlesen
        IF (.NOT. KSHELL(ION)) CYCLE
        DO
         READ (ION,'(A)',END=6) ZEILE2 
         IF (ZEILE2(:10) .EQ. 'K-SHELL ') THEN 
C         Falls es eine K-Shell-Zeile gibt, soll einmal der Header 
C         geschrieben werden
          IF (ION .EQ. 1) WRITE (99,'(A)') 
     >      '*KEYWORD--*****SY*I*<-K-SIGMA><-K-SEXPO>***<-K-EION->'
          WRITE (99,'(A)') ZEILE2(:IDX(ZEILE2)) !Zeile schreiben
         ENDIF 
         CYCLE !naechste Zeile
    6    EXIT !Ende der Datei erreicht
        END DO !Ende der Schleife fuer ein Ion
        CYCLE !naechstes Ion
      END DO !Ende der Schleife ueber alle Ionen
                
      
C     Alle Ionendateien wieder schliessen         
   9  DO ION=1, NION
       CLOSE (ION)
      END DO

   99 RETURN
   
   
C********  ERROR EXITS **************************************
  900 WRITE(0,*) 'ERROR: MORE LEVELS THAN DIMENSIONED'
      GOTO 990
      
  901 WRITE(0,*) 'ERROR: Cannot open file:'
      WRITE(0,*) FILENAME(:IDX(FILENAME))
      GOTO 990
      
  903 WRITE(0,*) 'ERROR: cannot read levelnumber'
      GOTO 990        
          
            
  990 WRITE(0,*) '*** The error occured in the following line:'
      WRITE(0,*) FEHLERZEILE(ION)(:IDX(FEHLERZEILE(ION)))
      WRITE(0,*) 'Element was '//ELEMENT
      STOP '*** ERROR DETECTED IN SUBROUTINE NEWDATOMION'
      
         
      END 
      SUBROUTINE SARGC(TEXT,N)

C**   Diese Funktion ist aequivalent zur entsprechenden C-Funktion.
C**   Sie ermittelt die Anzahl der Argumente in einem String.
C**   Das Parsen des Strings uebernimmt die Routine sargp.
C**   Dort sind auch die syntaktischen Regeln beschrieben.

      CHARACTER*(*) TEXT
      INTEGER AS,AE,N,I

      I=0 				! Do not look for any Argument
      CALL SARGP(TEXT,N,I,AS,AE)

      RETURN

      END
      SUBROUTINE SARGP(TEXT,N,I,AS,AE)

C**	Die Subroutine sargp zerlegt einen String.
C**	Ermittelt werden:
C**	n : die Anzahl der Argumente und, 
C**     	falls i in [1..n],
C**	as und ae : Start- und Endindex des i-ten Arguments.
C**	Geparsed wird nach folgenden Regel:
C**	Leerzeichen werden grunsaetzlich nicht beachtet
C**	(Ausnahmen siehe unten).
C**	Argumente werden durch Leerzeichen oder Komma oder '=' oder ':'
C**     getrennt.
C**	Wird ein Argument durch Leerzeichen und Komma getrennt,
C**     so gilt dies als eine Trennung.
C**	Folgen zwei Kommata ohne Argument, also hoechstens durch
C**	Leerzeichen getrennt, gilt dies als Leerargument.
C**	Zeichen zwischen zwei doppelten Anfuehrungszeichen gelten
C**	als ein Argument, auch wenn Leerzeichen oder Kommata
C**	enthalten sind. In diesem Fall werden die Anfuehrungszeichen
C**	als nicht zum Argument gehoerig betrachtet.
C**	Interpretaion der Rueckgabewerte:
C**	n: immer die Anzahl der Argumente
C**	as=-1 -> Es wurde kein i-tes Argument gefunden (i nicht in [1..n])
C**	as=0  -> Das i-te Argument war leer (z.B. ",,")
C**	sonst : text(as:ae) = i-tes Argument

        CHARACTER*(*) TEXT
        INTEGER N                  ! Out: Anzahl der Argumente
        INTEGER I                  ! In : gesuchtes Argument
        INTEGER AS,AE              ! Out: erste u. letzte Position des 
                                   !         Arguments im String

	INTEGER TL,TI
        INTEGER STATE

	TL=LEN(TEXT)

	state=0
C	state=0 -> kein Argument aktiv
C	state=1 -> Argumentende gefunden, naechstes Komma
C			ist   k e i n   Leerargument
C	state=2 -> normales Argument aktiv
C	state=3 -> Argument in '"' aktiv


        N=0
        AS=-1
	AE=TL

	DO 100 TI=1,TL

C** Go here looking for next start of an argument
        IF (STATE .EQ. 0) THEN
	   IF (TEXT(TI:TI) .EQ. ' ') GOTO 100
           IF ((TEXT(TI:TI) .EQ. ',')
     $     .OR.(TEXT(TI:TI) .EQ. '=')
     $     .OR.(TEXT(TI:TI) .EQ. ':')) THEN
                            N=N+1
                            IF (N .EQ. I) THEN
                                AS=0
                            ENDIF
                            GOTO 100
           ENDIF
           IF (TEXT(TI:TI) .EQ. '"') THEN
                           N=N+1
			   IF (N .EQ. I) AS=TI+1
			   STATE=3
			   GOTO 100
           ENDIF
	   STATE=2
	   N=N+1
           IF (N .EQ. I) AS=TI
	   GOTO 100
        ELSEIF (STATE .EQ. 1) THEN
	   IF (TEXT(TI:TI) .EQ. ' ') GOTO 100
           IF ((TEXT(TI:TI) .EQ. ',')
     $     .OR.(TEXT(TI:TI) .EQ. '=')
     $     .OR.(TEXT(TI:TI) .EQ. ':')) THEN
			   STATE=0
			   GOTO 100
           ENDIF
           IF (TEXT(TI:TI) .EQ. '"') THEN
                           N=N+1
                           IF (N .EQ. I) AS=TI+1
			   STATE=3
			   GOTO 100
           ENDIF
	   STATE=2
	   N=N+1
           IF (N .EQ. I) AS=TI
	   GOTO 100
        ELSEIF (STATE .EQ. 2) THEN
           IF (TEXT(TI:TI) .EQ. ' ') THEN
                          STATE=1
			  IF (N .EQ. I) AE=TI
			  GOTO 100
	   ENDIF
	   IF ((TEXT(TI:TI) .EQ. ',')
     $     .OR.(TEXT(TI:TI) .EQ. '=')
     $     .OR.(TEXT(TI:TI) .EQ. ':')) THEN
			  STATE=0
                          IF (N .EQ. I) AE=TI-1
			  GOTO 100
           ENDIF
           GOTO 100
        ELSE 
C***  ! IF (STATE .EQ. 3)
           IF (TEXT(TI:TI) .EQ. '"') THEN
                          STATE=1
			  IF (N .EQ. I) AE=TI-1
			  GOTO 100
           ENDIF
        ENDIF

100	CONTINUE

	RETURN
        END
      SUBROUTINE SARGV(TEXT,I,ARGTEXT)

C**   Diese Funktion ist fast aequivalent zur entsprechenden C-Funktion.
C**   Sie ermittelt das i-te Argument in einem String.
C**   Die Regeln zur Argumenttrennung sind in sargp beschrieben.
C**   Ist das i-te Argument nicht vorhanden, wird argtext nicht veraendert.

      CHARACTER*(*) TEXT,ARGTEXT
      INTEGER I
      INTEGER N,AS,AE

      CALL SARGP(TEXT,N,I,AS,AE)

      IF (AS .EQ. -1) GOTO 10
      IF (AS .EQ. 0) THEN
             ARGTEXT=' '
      ELSE
             ARGTEXT=TEXT(AS:AE)
      ENDIF

10    RETURN

      END

      PROGRAM MAINcoli 
C***  Provide Link data for possible use in the programm
      CHARACTER LINK_DATE*30, LINK_USER*10, LINK_HOST*60
      COMMON / COM_LINKINFO / LINK_DATE, LINK_USER, LINK_HOST
      LINK_DATE = ''
      LINK_USER = ''
      LINK_HOST = ''
                               
      CALL coli 
      END

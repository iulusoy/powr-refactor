      PROGRAM MAINnewformal_cards 
C***  Provide Link data for possible use in the programm
      CHARACTER LINK_DATE*30, LINK_USER*10, LINK_HOST*60
      COMMON / COM_LINKINFO / LINK_DATE, LINK_USER, LINK_HOST
      LINK_DATE = 'Di 21. Nov 13:00:30 CET 2023'
      LINK_USER = 'inga'
      LINK_HOST = 'ssc-laptop01'
                               
      CALL newformal_cards 
      END

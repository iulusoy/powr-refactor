
# Makefile to create all PoWR exe-files (.exe and .exe.opt) from the .f files in the same directory
# H. Todt 09.10.2022 htodt@astro.physik.uni-potsdam.de

linkdate=$(shell date)
linkuser=$(shell whoami)
linkhost=$(shell hostname)

FC = ifort -i8 -r8
FFLAGS     = -assume byterecl -save -extend-source  -O1 -fpe0 -traceback -mcmodel medium -g -fpconstant -fp-model strict
FFLAGS_OPT = -assume byterecl -save -extend-source  -O2 -fpe0 -traceback -mcmodel medium -g -axCORE-AVX2 -ipo -unroll -fpconstant -fp-model strict

MKLPATH    = ${MKLROOT}/lib/intel64
MKLINCLUDE = ${MKLROOT}/include

LINKER_OPTIONS = -L${MKLPATH} -I${MKLINCLUDE} -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -mcmodel medium 
LINKER_DYNAMIC = -shared-intel 
# LINKER_STATIC  = -L${MKLPATH} -I${MKLINCLUDE} -Wl,--start-group ${MKLPATH}/libmkl_intel_lp64.a ${MKLPATH}/libmkl_sequential.a ${MKLPATH}/libmkl_core.a -Wl,--end-group -lpthread 

COLIOBJECTS     =  addhistentry.o  addopa.o  append_autolevels.o  backjc.o  bnue.o  calc_s.o  check_cont.o  check_lines.o  cldiffus.o  clloade.o  clock.o  clopene.o  closms.o  clsavee.o  clsavewc.o  cmfcoop.o  cmffeop.o  cmsstore.o  coli.o  colihist.o  colimo.o  colimop.o  coli_setzero.o  coliwm.o  coop.o  count.o  datom.o  dbnuedt.o  deccoli.o  difdtdr.o  drlevel.o  fecheck.o  fedat.o  findcharge.o  frequint.o  frequbak.o  frequbakion.o  frequnorm.o  gauntff.o  genwp1.o  gethistentry.o  horner.o  idx.o  install.o  inv.o  isrcheq.o  isrchfgt.o  jsymset.o  ksigma.o  liop.o  lintridiagsol.o  opaross.o  openexms.o  openms.o  owninv.o  photocs.o  photon3.o  plotalpha.o  plotanf.o  plotcon.o  plotcons.o  plotrtau1coli.o  polyfit.o  popmin_nulling.o  preline.o  prep_drlines.o  prepk.o  prep_ppp.o  pri_epsg.o  priopacoli.o  readms.o  remark.o  rmodcoli.o  rsort.o  sargc.o  sargp.o  sargv.o  seqlinecl.o  set_momzero.o  setup_ff.o  sfit.o  sfitplo.o  shortray.o  splinpo_fast.o  splinpox.o  stamp.o  storage.o  store_ff.o  tradfun.o  trbk.o  vmf.o  vdopdd_setup.o  wmodcoli.o  writms.o 

.PHONY : all all_exe all_opt clean

all : all_exe all_opt
all_exe : coli.exe
all_opt : coli.exe.opt

coli.exe : $(COLIOBJECTS) maincoli.for
	$(FC) -o $(@) $(FFLAGS) $(LINKER_OPTIONS) $(LINKER_DYNAMIC) $(^)

coli.exe.opt : maincoli.for $(COLIOBJECTS) 
	$(shell cat maincoli.for $(patsubst %.o, %.f,$(COLIOBJECTS)) >> all.for )
	$(shell mv all.for maincoli.for)
	$(FC) -o $(@) $(FFLAGS_OPT)  $(LINKER_OPTIONS) $(LINKER_DYNAMIC) maincoli.for


# for each program a file main${progam}.for is required in which ${progam}.f is called
%.for :
	@printf "      PROGRAM MAIN$(patsubst main%.for,%,$(@)) \n"                > main$(patsubst main%.for,%,$(@)).for
	@printf "C***  Provide Link data for possible use in the programm\n"      >> main$(patsubst main%.for,%,$(@)).for
	@printf "      CHARACTER LINK_DATE*30, LINK_USER*10, LINK_HOST*60\n"      >> main$(patsubst main%.for,%,$(@)).for
	@printf "      COMMON / COM_LINKINFO / LINK_DATE, LINK_USER, LINK_HOST\n" >> main$(patsubst main%.for,%,$(@)).for
	@printf "      LINK_DATE = '$(linkdate)'\n"                               >> main$(patsubst main%.for,%,$(@)).for
	@printf "      LINK_USER = '$(linkuser)'\n"                               >> main$(patsubst main%.for,%,$(@)).for
	@printf "      LINK_HOST = '$(linkhost)'\n"                               >> main$(patsubst main%.for,%,$(@)).for
	@printf "                               \n"                               >> main$(patsubst main%.for,%,$(@)).for
	@printf "      CALL $(patsubst main%.for,%,$(@)) \n"                      >> main$(patsubst main%.for,%,$(@)).for
	@printf "      END\n"                                                     >> main$(patsubst main%.for,%,$(@)).for


.f.o :
	$(FC) $(FFLAGS) -c $?

clean :
	rm -f *.o
	rm -f *.for
	rm -f *.exe
	rm -f *.exe.opt

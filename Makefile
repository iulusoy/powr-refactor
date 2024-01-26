# Makefile for building PoWR
# ISUlusoy, SSC, 11/23

# necessary directories
SRC_DIR = src
OBJ_DIR = build
BIN_DIR = powr/exe.dir
BIN_DIR_DEBUG = powr/exe_dev.dir
LIB_DIR = lib

# compiler and linking options
FC = ifort
FFLAGS = -i8 -r8 -I${LIB_DIR} -assume byterecl -save -extend-source -O3 -fpe0 -traceback -mcmodel medium -g -fpconstant -fp-model strict
MKLPATH    = ${MKLROOT}/lib/intel64
MKLINCLUDE = ${MKLROOT}/include
LINKER_OPTIONS = -L${MKLPATH} -I${MKLINCLUDE} -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -mcmodel medium
LINKER_DYNAMIC = -shared-intel

# relevant source files for the different programs
ADAPTERSRC = adapop.f adapter.f adatrans.f addhistentry.f append_autolevels.f chrinstr.f clock.f closms.f cmsstore.f count.f datom.f decadp.f fedat.f findcharge.f idx.f install.f isrcheq.f jsymset.f lengthms.f lipo.f openms.f priadp.f readms.f remark.f rmodadp.f sargc.f sargp.f sargrest.f sargv.f second.f splinpox.f stamp.f storage.f trbk.f writms.f mainadapter.f
ADAPTERSRCDIR = $(addprefix $(SRC_DIR)/,$(ADAPTERSRC))
ADAPTEROBJ = $(patsubst $(SRC_DIR)/%.f, $(OBJ_DIR)/%.o, $(ADAPTERSRCDIR))

COMOSRC = addhistentry.f append_autolevels.f bfcross.f bnue.f clock.f closms.f cmsstore.f como.f coop.f count.f datom.f dbnuedt.f decnot.f decomo.f difdtdr.f diffus.f fedat.f findcharge.f folr.f frequbak.f frequbakion.f gauntff.f gethistentry.f idx.f install.f invtri.f isrcheq.f jsymset.f ksigma.f momo.f opaross.f openms.f photocs.f photon3.f plohtot.f plotanf.f plotanfs.f plotcon.f plotcons.f plotopa.f plotrtau1.f popmin_nulling.f prihtot.f primint.f priopa.f readms.f remark.f remoco.f sargc.f sargp.f sargv.f second.f stamp.f storage.f tradfun.f trbk.f writms.f maincomo.f
COMOSRCDIR = $(addprefix $(SRC_DIR)/,$(COMOSRC))
COMOOBJ = $(patsubst $(SRC_DIR)/%.f, $(OBJ_DIR)/%.o, $(COMOSRCDIR))

COLISRC = addhistentry.f addopa.f append_autolevels.f backjc.f bnue.f calc_s.f check_cont.f check_lines.f cldiffus.f clloade.f clock.f clopene.f closms.f clsavee.f clsavewc.f cmfcoop.f cmffeop.f cmsstore.f coli.f colihist.f colimo.f colimop.f coli_setzero.f coliwm.f coop.f count.f datom.f dbnuedt.f deccoli.f difdtdr.f drlevel.f fecheck.f fedat.f findcharge.f folr.f frequint.f frequbak.f frequbakion.f frequnorm.f gauntff.f genwp1.f gethistentry.f horner.f idx.f install.f inv.f isamax.f isrcheq.f isrchfgt.f jsymset.f ksigma.f liop.f lintridiagsol.f opaross.f openexms.f openms.f owninv.f photocs.f photon3.f plotalpha.f plotanf.f plotcon.f plotcons.f plotrtau1coli.f polyfit.f popmin_nulling.f preline.f prep_drlines.f prepk.f prep_ppp.f pri_epsg.f priopacoli.f readms.f remark.f rmodcoli.f rsort.f sargc.f sargp.f sargv.f second.f seqlinecl.f set_momzero.f setup_ff.f sfit.f sfitplo.f shift.f shortray.f splinpo_fast.f splinpox.f stamp.f storage.f store_ff.f tradfun.f trbk.f vmf.f vdopdd_setup.f wmodcoli.f writms.f maincoli.f
COLISRCDIR = $(addprefix $(SRC_DIR)/,$(COLISRC))
COLIOBJ = $(patsubst $(SRC_DIR)/%.f, $(OBJ_DIR)/%.o, $(COLISRCDIR))

EXTRAPSRC = addhistentry.f aitken.f append_autolevels.f change.f clock.f closms.f cmsstore.f count.f datom.f decnot.f extrap.f fedat.f findcharge.f gauntff.f gethistentry.f idx.f inhibit.f isrcheq.f jsymset.f ng3.f ng4.f openms.f pricorr.f priex.f readms.f remark.f sargc.f sargp.f sargv.f second.f stamp.f storage.f writms.f mainextrap.f
EXTRAPSRCDIR = $(addprefix $(SRC_DIR)/,$(EXTRAPSRC))
EXTRAPOBJ = $(patsubst $(SRC_DIR)/%.f, $(OBJ_DIR)/%.o, $(EXTRAPSRCDIR))

FORMALSRC = append_autolevels.f backjc.f bandwidth.f bnue.f clock.f closms.f cmffeop.f cmfset_formal.f cmsstore.f convolgauss_flex.f convolopafe.f coop.f copy_secondmodel.f count.f datom.f dbnuedt.f decform.f diffus.f drtrans.f elimin.f equal.f fecheck.f fedat.f fierfc.f filterfunctions.f findcharge.f findind.f findldr.f folr.f formal.f formcmf.f formosa.f gauntff.f genw0.f gmalu.f horner.f idx.f insert_line.f install.f inv.f invtri.f isamax.f isrcheq.f isrchfge.f isrchfgt.f isrchfle.f isrchflt.f jsymset.f kholtsmark.f ksigma.f limbdark_output.f limbdark_prep.f linstark.f liop.f lipo.f ltepop.f macroclump.f manipop.f mdmv.f mdv.f merge_rgrid.f moment0.f moment1.f moment2.f msub.f multiple.f multspli.f mvmd.f mvv.f newpol2.f nowind.f obsfram.f openms.f owninv.f phiholtsmark.f photocs.f photon3.f plotanf.f plotanfs.f plotcon.f plotcons.f plot_secondmodel_grid.f plotvdop.f plot_windrot_grid.f polyfit.f popmin_nulling.f preform.f prepmacroclump.f prepray.f pridwl.f priopal.f pri_par.f pripro.f quadstark.f read_h_starkdata.f read_linecard_parameters.f readms.f remark.f rescale_secmod.f rotation_prep.f sargc.f sargp.f sargrest.f sargv.f sdot.f second.f secondmodel_define.f secondmodel_prep.f set_pop_zero.f setup.f sfit.f sfitplo.f shift.f sofbet.f splinpo.f splinpo_fast.f splinpox.f stamp.f starkbroad.f starkdamp_hei.f starkheii.f starkheiiprep.f starkhi.f stark_hi_lemke.f starkholtsmark.f starkprof.f starkvoigt.f storage.f taucmu.f tradfun.f tradwl.f transform_rgrid.f traplo.f trbk.f vadd.f vcse1f.f vdop_struct.f vmalv.f vmf.f voigth.f writms.f zoneint.f mainformal.f
FORMALSRCDIR = $(addprefix $(SRC_DIR)/,$(FORMALSRC))
FORMALOBJ = $(patsubst $(SRC_DIR)/%.f, $(OBJ_DIR)/%.o, $(FORMALSRCDIR))

MODIFYSRC = addhistentry.f append_autolevels.f change.f clock.f closms.f cmsstore.f count.f datom.f fedat.f findcharge.f gethistentry.f idx.f inhibit.f intepo.f isrcheq.f jsymset.f modify.f openms.f primo.f readms.f remark.f sargc.f sargp.f sargv.f second.f stamp.f storage.f writms.f mainmodify.f
MODIFYSRCDIR = $(addprefix $(SRC_DIR)/,$(MODIFYSRC))
MODIFYOBJ = $(patsubst $(SRC_DIR)/%.f, $(OBJ_DIR)/%.o, $(MODIFYSRCDIR))

MSINFOSRC = idx.f msinfo.f storage.f mainmsinfo.f
MSINFOSRCDIR = $(addprefix $(SRC_DIR)/,$(MSINFOSRC))
MSINFOOBJ = $(patsubst $(SRC_DIR)/%.f, $(OBJ_DIR)/%.o, $(MSINFOSRCDIR))

NEWDATOMSRC = findelement.f idx.f jsymset.f newdatom.f newdatomion.f sargc.f sargp.f sargv.f mainnewdatom.f
NEWDATOMSRCDIR = $(addprefix $(SRC_DIR)/,$(NEWDATOMSRC))
NEWDATOMOBJ = $(patsubst $(SRC_DIR)/%.f, $(OBJ_DIR)/%.o, $(NEWDATOMSRCDIR))

NEWFORMAL_CARDSSRC = append_autolevels.f clock.f closms.f cmsstore.f count.f datom.f fedat.f findcharge.f idx.f install.f isrcheq.f jsymset.f newformal_cards.f openms.f readms.f remark.f sargc.f sargp.f sargrest.f sargv.f storage.f mainnewformal_cards.f
NEWFORMAL_CARDSSRCDIR = $(addprefix $(SRC_DIR)/,$(NEWFORMAL_CARDSSRC))
NEWFORMAL_CARDSOBJ = $(patsubst $(SRC_DIR)/%.f, $(OBJ_DIR)/%.o, $(NEWFORMAL_CARDSSRCDIR))

NJNSRC = addhistentry.f append_autolevels.f closms.f cmsstore.f idx.f jsymset.f njn.f openms.f readms.f storage.f writms.f mainnjn.f
NJNSRCDIR = $(addprefix $(SRC_DIR)/,$(NJNSRC))
NJNOBJ = $(patsubst $(SRC_DIR)/%.f, $(OBJ_DIR)/%.o, $(NJNSRCDIR))

STEALSRC = addhistentry.f addopa.f adjgamma.f aitken.f append_autolevels.f backuppopnum.f bfcross.f bnue.f brmv.f brnorm2.f brtpup.f brvdivs.f brvm.f brvvdy.f calcgammaradmean.f calcmassfromgeff.f calcwc.f cbbfe.f cbbh.f cbbhe.f cbbmore.f cbbn.f ccore.f change.f clock.f closms.f clump_struct.f cmsstore.f cofreq.f colli.f coma.f coop.f coopfrq.f count.f datom.f dbnuedt.f dcoop.f decnot.f decste.f decvelpar.f delpla.f deltagr.f deltagrthin.f deriv.f dliop.f dmopen.f drdat.f ensuretaumax.f erf.f expint1exp.f extrap.f fecheck.f fedat.f feop_steal.f filterfunctions.f findcharge.f flag_zerorates.f flgrid.f funsch.f gauntff.f geomesh.f gethistentry.f gradiff.f hysthdruku.f hystruku.f idx.f inhibit.f initfcorr.f initvel.f initvelbetapar.f install.f interpolatepopnum.f interpolatetemp.f inv.f isamax.f ismax.f isrcheq.f isrchfge.f isrchfgt.f isrchfle.f isrchflt.f isrchne.f jlderiv.f jsymset.f ksigma.f lcore.f lengthms.f linpop.f linsol.f linsol_split.f liop.f lipo.f load_ff.f loadwc.f ltepop.f mgoetz.f mlanger.f nextjob.f ng3.f ng4.f nltepop.f opaross.f openexms.f openms.f overlap.f owninv.f pgrid.f photocs.f photon3.f plocc.f plotacc.f plotaccelem.f plotalpha.f plotanf.f plotanfs.f plotapp.f plotcon.f plotcons.f plotdep.f plotflu.f plotfgadapter.f plotfgstrat.f plotforcemult.f plotgamma.f plothsum.f plotjline.f plotjnue.f plotpop.f plotsigmafe.f plott.f plotunlu.f plotv.f plotvgrad.f pop_renorm.f popsmall.f popzero.f preline.f prep_drlines.f prepkubat.f preplotapp.f pri1rat.f pricc.f pricolr.f pricomp.f pricorr.f pridat.f priex.f priexpo.f priflux.f prigahist.f prihist.f prijost.f prilc.f primat.f primod.f printmodelsummary.f pripop.f prirat.f pritau.f priunlu.f radio.f radnet.f readms.f redcor.f regridoldpop.f regula.f remark.f remarkf.f remost.f rgrid.f rsort.f sargc.f sargp.f sargv.f scaledm.f sdot.f second.f seqlinecl.f setxjc.f setxjfine.f setxjl.f setxjlcf.f shift.f shiftreal.f shiftstring.f smach.f splinpo.f splinpox.f stamp.f steal.f sthist.f storage.f tauscal.f tcolor.f tdiffus.f tempcorr.f tempcorr_expdamp.f tempcorr_fluxerr.f tradfun.f trbk.f vdopdd_setup.f velobeta.f velthin.f vmf.f vsub.f writms.f wrvel.f xextinc.f xrudi.f zanstra.f mainsteal.f
STEALSRCDIR = $(addprefix $(SRC_DIR)/,$(STEALSRC))
STEALOBJ = $(patsubst $(SRC_DIR)/%.f, $(OBJ_DIR)/%.o, $(STEALSRCDIR))

WRCONTSRC = addhistentry.f append_autolevels.f bnue.f clock.f closms.f cmsstore.f coop.f count.f datom.f dbnuedt.f decon.f delpla.f difdtdr.f diffus.f elimin.f equal.f fedat.f filterfunctions.f findcharge.f gauntff.f gethistentry.f horner.f idx.f install.f inv.f isamax.f isrcheq.f jsymset.f ksigma.f lipo.f mdmv.f mdv.f moment0.f moment1.f moment2.f msub.f mvmd.f mvv.f opaross.f openms.f openmsr.f owninv.f photocs.f photon3.f plotanf.f plotanfs.f plotcon.f plotcons.f polyfit.f popmin_nulling.f pricolr.f priflux.f priint.f prijost.f priopa.f readms.f regula.f remark.f rmodcon.f sargc.f sargp.f sargv.f second.f setup.f sfit.f sfitplo.f splinpo.f splinpox.f stamp.f storage.f tcolor.f tradfun.f trbk.f tremain.f vadd.f vmf.f wrcont.f writms.f xextinc.f xrudi.f zanstra.f mainwrcont.f
WRCONTSRCDIR = $(addprefix $(SRC_DIR)/,$(WRCONTSRC))
WRCONTOBJ = $(patsubst $(SRC_DIR)/%.f, $(OBJ_DIR)/%.o, $(WRCONTSRCDIR))

WRSTARTSRC = addhistentry.f append_autolevels.f bfcross.f bnue.f clock.f closms.f clump_struct.f cmsstore.f coop.f coopfrq.f count.f datom.f dbnuedt.f dec2beta.f decfreq.f decstar.f decvelpar.f deltagr.f deltagrthin.f dtdr.f fedat.f fgrid.f findcharge.f gauntff.f geomesh.f gradiff.f grey.f hysthdruku.f hystruku.f idx.f initvel.f initvelbetapar.f install.f isrcheq.f isrchfgt.f isrchne.f jstart.f jsymset.f ksigma.f lipo.f loadoldj.f lowercase.f ltepop.f mgoetz.f mlanger.f my_clock.f my_date.f opagrey.f opaross.f openms.f pgrid.f photocs.f photon3.f plotanf.f plotanfs.f plotcon.f plotcons.f plotv.f plotvgrad.f prep_drlines.f prep_gammarad.f pricomp.f primod.f priparam.f prixdat.f readms.f readoldt.f regula.f remark.f rgrid.f ruku.f sargc.f sargp.f sargv.f second.f sequin.f shiftreal.f shiftstring.f splinpo.f splinpox.f stamp.f storage.f tabread.f tauscal.f tradfun.f trbk.f velobeta.f velthin.f vturb_setup.f writms.f wrstart.f wrvel.f xrudi.f mainwrstart.f
WRSTARTSRCDIR = $(addprefix $(SRC_DIR)/,$(WRSTARTSRC))
WRSTARTOBJ = $(patsubst $(SRC_DIR)/%.f, $(OBJ_DIR)/%.o, $(WRSTARTSRCDIR))

# the different targets
all: adapter como coli extrap formal modify msinfo newdatom newformal_cards njn steal wrcont wrstart

adapter: adapter.exe
coli: coli.exe
como: como.exe
extrap: extrap.exe
formal: formal.exe
modify: modify.exe
msinfo: msinfo.exe
newdatom: newdatom.exe
newformal_cards: newformal_cards.exe
njn: njn.exe
steal: steal.exe
wrcont: wrcont.exe
wrstart: wrstart.exe

# debug options and rules
# debug: FFLAGS = -i8 -r8 -I${LIB_DIR} -assume byterecl -save -extend-source -O0 -fpe0 -traceback -mcmodel medium -g -fpconstant -fp-model strict -warn all -check all -traceback -fp-stack-check
debug: FFLAGS = -i8 -r8 -I${LIB_DIR} -assume byterecl -save -extend-source -O0 -fpe0 -traceback -mcmodel medium -g -fpconstant -fp-model strict -check all -traceback -fp-stack-check
debug: coli steal wrcont como
debug: BIN_DIR = $(BIN_DIR_DEBUG)

# debug_all: FFLAGS = -i8 -r8 -I${LIB_DIR} -assume byterecl -save -extend-source -O0 -fpe0 -traceback -mcmodel medium -g -fpconstant -fp-model strict
debug_all: FFLAGS = -i8 -r8 -I${LIB_DIR} -assume byterecl -save -extend-source -O0 -fpe0 -traceback -mcmodel medium -g -fpconstant -fp-model strict
debug_all: BIN_DIR = $(BIN_DIR_DEBUG)
debug_all: adapter como coli extrap formal modify msinfo newdatom newformal_cards njn steal wrcont wrstart

adapter.exe: $(ADAPTEROBJ)
	$(FC) $(FFLAGS) $(LINKER_OPTIONS) $(LINKER_DYNAMIC) -o $(BIN_DIR)/$@ $^

como.exe: $(COMOOBJ)
	$(FC) $(FFLAGS) $(LINKER_OPTIONS) $(LINKER_DYNAMIC) -o $(BIN_DIR)/$@ $^

coli.exe: $(COLIOBJ)
	$(FC) $(FFLAGS) $(LINKER_OPTIONS) $(LINKER_DYNAMIC) -o $(BIN_DIR)/$@ $^

extrap.exe: $(EXTRAPOBJ)
	$(FC) $(FFLAGS) $(LINKER_OPTIONS) $(LINKER_DYNAMIC) -o $(BIN_DIR)/$@ $^

formal.exe: $(FORMALOBJ)
	$(FC) $(FFLAGS) $(LINKER_OPTIONS) $(LINKER_DYNAMIC) -o $(BIN_DIR)/$@ $^

modify.exe: $(MODIFYOBJ)
	$(FC) $(FFLAGS) $(LINKER_OPTIONS) $(LINKER_DYNAMIC) -o $(BIN_DIR)/$@ $^

msinfo.exe: $(MSINFOOBJ)
	$(FC) $(FFLAGS) $(LINKER_OPTIONS) $(LINKER_DYNAMIC) -o $(BIN_DIR)/$@ $^

newdatom.exe: $(NEWDATOMOBJ)
	$(FC) $(FFLAGS) $(LINKER_OPTIONS) $(LINKER_DYNAMIC) -o $(BIN_DIR)/$@ $^

newformal_cards.exe: $(NEWFORMAL_CARDSOBJ)
	$(FC) $(FFLAGS) $(LINKER_OPTIONS) $(LINKER_DYNAMIC) -o $(BIN_DIR)/$@ $^

njn.exe: $(NJNOBJ)
	$(FC) $(FFLAGS) $(LINKER_OPTIONS) $(LINKER_DYNAMIC) -o $(BIN_DIR)/$@ $^

steal.exe: $(STEALOBJ)
	$(FC) $(FFLAGS) $(LINKER_OPTIONS) $(LINKER_DYNAMIC) -o $(BIN_DIR)/$@ $^

wrcont.exe: $(WRCONTOBJ)
	$(FC) $(FFLAGS) $(LINKER_OPTIONS) $(LINKER_DYNAMIC) -o $(BIN_DIR)/$@ $^

wrstart.exe: $(WRSTARTOBJ)
	$(FC) $(FFLAGS) $(LINKER_OPTIONS) $(LINKER_DYNAMIC) -o $(BIN_DIR)/$@ $^

# rules to compile
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f
	$(FC) $(FFLAGS) -c $< -o $@

print_info:
	$(info COLISRC: $(COLISRC))
	$(info COLIOBJ: $(COLIOBJ))

clean:
	rm -f build/*.o $(BIN_DIR)/*.exe $(BIN_DIR_DEBUG)/*.exe

clean_build:
	rm -f build/*.o
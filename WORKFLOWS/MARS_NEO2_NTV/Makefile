NEO2PATH = $(HOME)/PHD/INTEGRATED_DEVELOPMENT_REPO/NEO-2
EXECUTABLEPATH = $(NEO2PATH)/NEO-2-QL/Build/

MATLABPATH = $(NEO2PATH)/OctaveScripts/

MARS_DIR = '/proj/plasma/DATA/DEMO/teams/MARSQ_OUTPUTS/DATA_equil_qmod/MARSQ_OUTPUTS_100kAt_dBkinetic_NTVkinetic_NEO2profs'
OUTPUT_DIR = './input_files_for_generate_neo2_profile/'

PATH_TO_INPUTFILES = $(OUTPUT_DIR)
NUMBER_OF_SURFACES = 10
BOUNDS = [0.693,0.695] # resonant surface at stor = 0.694


AXIPATH = /temp/ert/CONDOR/Neo2/DEMO/fromefit_poscur_negflux.bc
PERTPATH = /proj/plasma/DATA/DEMO/BOOZER/frommars-n_minus1_v4.bc

USER_MAIL = georg.grassler@tugraz.at
REQUESTED_RAM = 29*1024

############################################################################################################

all: runs

runs: run_folders submit submit_failed.template
	condor_submit submit

run_folders: neo_2.x multi_spec_demo.in neo.in neo2.in axi.bc pert.bc
	./neo_2.x

neo_2.x:
	cp $(EXECUTABLEPATH)/neo_2.x neo_2.x

# Deprecated Matlab version
# multi_spec_demo.in: input_files_for_generate_neo2_profile
# 	matlab -nodisplay -r "addpath('$(MATLABPATH)'); run_generate_neo2_profile($(NUMBER_OF_SURFACES),$(BOUNDS)); exit;"

multi_spec_demo.in: input_files_for_generate_neo2_profile
	python3 -c "\
	from neo2_mars import run_generate_neo2_profile; \
	run_generate_neo2_profile($(PATH_TO_INPUTFILES), $(NUMBER_OF_SURFACES), $(BOUNDS)); \
	"

input_files_for_generate_neo2_profile:
	python3 -c "\
	from neo2_mars import write_input_for_generate_neo2_profile_from_mars; \
	write_input_for_generate_neo2_profile_from_mars($(MARS_DIR), $(OUTPUT_DIR)); \
	"

neo.in:
	cp $(NEO2PATH)/DOC/neo.in.ql-full neo.in
	sed -i "s/example_bg_field\.bc/axi.bc/" neo.in

neo2.in:
	cp $(NEO2PATH)/DOC/neo2.in.ql-full neo2.in
	sed -i "s/\(^\s*in_file_pert\s*=\s*'\)[^']*'/\1pert.bc'/" neo2.in
	sed -i 's/\(^\s*fname_multispec_in\s*=\s*"\)[^"]*"/\1multi_spec_demo.in"/' neo2.in
	sed -i 's/isw_multispecies_init\s*=\s*[^ ]*/isw_multispecies_init = 1/' neo2.in

axi.bc:
	ln -s $(AXIPATH) axi.bc

pert.bc:
	ln -s $(PERTPATH) pert.bc


submit:
	cp $(NEO2PATH)/DOC/condor_submit_example submit
	sed -i '19 i environment = "XDG_CACHE_HOME=/tmp OMP_NUM_THREADS=1"' submit
	sed -i 's/notify_user\s*=\s*[^ ]*/notify_user = $(USER_MAIL)/' submit
	sed -i 's/request_memory\s*=\s*[^ ]*/request_memory = $(REQUESTED_RAM)/' submit

submit_failed.template:
	cp $(NEO2PATH)/DOC/condor_submit_failed.template submit_failed.template
	sed -i 's/notify_user\s*=\s*[^ ]*/notify_user = $(USER_MAIL)/' submit_failed.template
	

clean:
	-rm neo.in
	-rm neo2.in
	-rm submit
	-rm submit_failed.template
	-rm *.bc
	-rm *.dat
	-rm -r input_files_for_generate_neo2_profile
	-rm multi_spec_demo.in
	-rm *.h5
	-rm -r es_*
	-rm err.*
	-rm *.log
#	rm .in *.bc

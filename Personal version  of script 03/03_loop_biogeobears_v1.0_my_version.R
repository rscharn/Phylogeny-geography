#This script runs biogeobears on all tree/geog files in a folder (provided that they are in the same format of the script I made to generate them)
#######################INPUTS######################################

#set working directory (where all imput data is present, including shapefiles)
setwd("C:\\Users\\Ruud\\Desktop\\R working Directory\\scripts fernanda\\LASTEST VERSIONS\\Ferns")
require(optimx)         # You need to have some version of optimx available
                        # as it is a BioGeoBEARS dependency; however, if you
                        # don't want to use optimx, and use optim() (from R core) 
                        # you can set:
                        # BioGeoBEARS_run_object$use_optimx = FALSE
                        # ...everything should work either way -- NJM 2014-01-08
require(FD)       # for FD::maxent() (make sure this is up-to-date)
#library(snow)     # (if you want to use multicore functionality; some systems/R versions prefer library(parallel), try either)
require(parallel)

require(BioGeoBEARS)
#source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_basics_v1.R")
#source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_classes_v1.R")
#source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_models_v1.R")
#source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_plots_v1.R")
#source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_readwrite_v1.R")
#source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_simulate_v1.R")
#source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_stratified_v1.R")
#source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_univ_model_v1.R")
#source("http://phylo.wikidot.com/local--files/biogeobears/calc_loglike_sp_v01.R")
source ("C:\\Users\\Ruud\\Desktop\\R working Directory\\Fossil analyses\\Biogeobears Important dependancys\\BioGeoBEARS_basics_v1.R")
source ("C:\\Users\\Ruud\\Desktop\\R working Directory\\Fossil analyses\\Biogeobears Important dependancys\\BioGeoBEARS_classes_v1.R")
source ("C:\\Users\\Ruud\\Desktop\\R working Directory\\Fossil analyses\\Biogeobears Important dependancys\\BioGeoBEARS_models_v1.R")
source ("C:\\Users\\Ruud\\Desktop\\R working Directory\\Fossil analyses\\Biogeobears Important dependancys\\BioGeoBEARS_plots_v1.R")
source ("C:\\Users\\Ruud\\Desktop\\R working Directory\\Fossil analyses\\Biogeobears Important dependancys\\BioGeoBEARS_readwrite_v1.R")
source ("C:\\Users\\Ruud\\Desktop\\R working Directory\\Fossil analyses\\Biogeobears Important dependancys\\BioGeoBEARS_simulate_v1.R")
source ("C:\\Users\\Ruud\\Desktop\\R working Directory\\Fossil analyses\\Biogeobears Important dependancys\\BioGeoBEARS_stratified_v1.R")
source ("C:\\Users\\Ruud\\Desktop\\R working Directory\\Fossil analyses\\Biogeobears Important dependancys\\BioGeoBEARS_univ_model_v1.R")
source ("C:\\Users\\Ruud\\Desktop\\R working Directory\\Fossil analyses\\Biogeobears Important dependancys\\calc_loglike_sp_v01.R")

pdf = TRUE #do you want pdf's made visualizing the reconstruction for each tree
###################################################################
biogeobears.loop <-function(wd=getwd(),max_range_size = 5,pdf = TRUE){

calc_loglike_sp = compiler::cmpfun(calc_loglike_sp_prebyte)    # crucial to fix bug in uppass calculations
calc_independent_likelihoods_on_each_branch = compiler::cmpfun(calc_independent_likelihoods_on_each_branch_prebyte)

number=length(grep("clade_\\d+_area_\\d\\_updated\\.tre",dir(),perl =TRUE))

for (i in 1:number){
	#i = 33
	a = "clade_"  
	c = "_area_1.geog"
	d = "_area_1_updated.tre"
	e = "area_1_DEC"
	f = ".Rdata"
	g = ".pdf"
		
	#create 50 geogs
	treeobjectname=paste(a,i,d,sep="")
	geogobjectname=paste(a,i,c,sep="")
	rdataobject=paste(a,i,e,f,sep="")
	pdfname=paste(a,i,e,g,sep="")
	
trfn = treeobjectname
geogfn = geogobjectname

moref(trfn)

# Look at your phylogeny:
tr = read.tree(trfn)
tr

# Look at the raw geography text file:
moref(geogfn)

# Look at your geographic range data:
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
tipranges
tail(tipranges@df)


# Set the maximum number of areas any species may occupy; this cannot be larger 
# than the number of areas you set up, but it can be smaller.


# Set up empty tables to hold the statistical results
restable = NULL
teststable = NULL

#######################################################
# Run DEC
#######################################################

BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$use_optimx = TRUE
BioGeoBEARS_run_object$speedup=TRUE        # seems to work OK
BioGeoBEARS_run_object$calc_ancprobs=TRUE        # get ancestral states from optim run

# Set up a time-stratified analysis 
# (un-comment to use; see example files in extdata_dir, 
#  and BioGeoBEARS google group posts for further hints)
#BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
#BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
#BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"

# Input the maximum range size
BioGeoBEARS_run_object$max_range_size = max_range_size

# Multicore processing if desired
BioGeoBEARS_run_object$num_cores_to_use=1
# (use more cores to speed it up; this requires
# library(parallel), which is default on Macs in
# R 3.0+, but apparently still has to be typed
# on Windows machines. Note: apparently parallel
# works on Mac command-line R, but not R.app.
# BioGeoBEARS checks for this and resets to 1
# core with R.app)

# Sparse matrix exponentiation is an option for huge numbers of ranges/states (600+)
# I have experimented with sparse matrix exponentiation in EXPOKIT/rexpokit,
# but the results are imprecise and so I haven't explored it further.
# In a Bayesian analysis, it might work OK, but the ML point estimates are
# not identical.
# Also, I have not implemented all functions to work with force_sparse=TRUE.
# Volunteers are welcome to work on it!!
BioGeoBEARS_run_object$force_sparse=FALSE

# Give BioGeoBEARS the location of the geography text file
BioGeoBEARS_run_object$geogfn = geogfn

# Give BioGeoBEARS the location of the phylogeny Newick file
BioGeoBEARS_run_object$trfn = trfn

# Read in the input files. Read the error messages if you get them!
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
# Divide the tree up by strata (stratified analysis)
#BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE
BioGeoBEARS_run_object$master_table

# Set up DEC model
# (nothing to do; defaults)

# Look at the BioGeoBEARS_run_object; it's just a list of settings etc.
BioGeoBEARS_run_object

# This contains the model object
BioGeoBEARS_run_object$BioGeoBEARS_model_object

# This table contains the parameters of the model 
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table

# Run this to check inputs. Read the error messages if you get them!
check_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Adding useAmbiguities (so that "?" can be used -- 
# note that this can lead to widespread ranges)
BioGeoBEARS_run_object$useAmbiguities = TRUE

# For a slow analysis, run once, then set runslow=FALSE to just 
# load the saved result.
runslow = TRUE
resfn = rdataobject
if (runslow)
    {
    res = bears_optim_run(BioGeoBEARS_run_object)
    res    

    save(res, file=resfn)
    resDEC = res
    } else {
    # Loads to "res"
    load(resfn)
    resDEC = res
    }

if(pdf==TRUE){
#######################################################
# PDF plots
#######################################################
pdffn = pdfname
pdf(pdffn, paper="a4")

#######################################################
# Plot ancestral states - DEC
#######################################################
analysis_titletxt ="BioGeoBEARS DEC"

# Setup
results_object = resDEC
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

dev.off()


#######################################################
# Statistics
#######################################################
# We have to extract the log-likelihood differently, depending on the 
# version of optim/optimx
if (BioGeoBEARS_run_object$use_optimx == TRUE)
    {
    # Using optimx() results
    if (packageVersion("optimx") < 2013)
        {
        # optimx 2012
        LnL_2 = as.numeric(resDEC$optim_result$fvalues)
        } else {
        # optimx 2013
        LnL_2 = as.numeric(resDEC$optim_result$value)
        } # end optimx 2012 vs. 2013
    } else {
    # Using optim() results
    LnL_2 = as.numeric(resDEC$optim_result$value)
    } # end optim vs. optimx

numparams1 = 3

res    # DEC, null model for Likelihood Ratio Test (LRT)
}
}
}

suppressWarnings(biogeobears.loop(wd=wd,pdf=pdf))
""" 
Module to set up run time parameters for Clawpack -- AMRClaw code.

The values set in the function setrun are then written out to data files
that will be read in by the Fortran code.

=====================================================================
Things that might change from one event to another:

    tfinal:  When to stop computation
    nout:  How many outputs frames 

    topofiles: Need to cover source region at 1 minute resolution and
               relevant part of ocean to at least 4 minute resolution.
    regions: Where to enforce refinement as wave propagates.
  
    amr_levels_max : Number of refinement levels to use (6 for inundation of CC)

    t_shelf, t_harbor: Times used in turning on gauges and regions
    tstart_max:   When to start tracking maxima on the fixed grid, 
                  Now uses max(100., t_harbor) by default
=====================================================================
    
""" 

import os
import numpy as np


t_shelf = 14.5*3600.   # time approaching continental slope
t_harbor = 15.*3600.  # time approaching harbor   # also used for fgmax below


#------------------------------
def setrun(claw_pkg='geoclaw'):
#------------------------------
    
    """ 
    Define the parameters used for running Clawpack.

    INPUT:
        claw_pkg expected to be "geoclaw" for this setrun.

    OUTPUT:
        rundata - object of class ClawRunData 
    
    """ 
    
    from clawpack.clawutil import data 
    
    
    assert claw_pkg.lower() == 'geoclaw',  "Expected claw_pkg = 'geoclaw'"

    num_dim = 2
    rundata = data.ClawRunData(claw_pkg, num_dim)

    #------------------------------------------------------------------
    # Problem-specific parameters to be written to setprob.data:
    #------------------------------------------------------------------



    #------------------------------------------------------------------
    # GeoClaw specific parameters:
    #------------------------------------------------------------------
    rundata = setgeo(rundata)
    
    #------------------------------------------------------------------
    # Standard Clawpack parameters to be written to claw.data:
    #------------------------------------------------------------------

    clawdata = rundata.clawdata  # initialized when rundata instantiated


    # Set single grid parameters first.
    # See below for AMR parameters.


    # ---------------
    # Spatial domain:
    # ---------------

    # Number of space dimensions:
    clawdata.num_dim = num_dim
    
    # Lower and upper edge of computational domain:
    clawdata.lower[0] = 215
    clawdata.upper[0] = 218
    clawdata.lower[1] = -19          # ylower
    clawdata.upper[1] = -16.0         # yupper
    
    # Number of grid cells:
    clawdata.num_cells[0] = 45     # mx
    clawdata.num_cells[1] = 45     # my

    # ---------------
    # Size of system:
    # ---------------

    # Number of equations in the system:
    clawdata.num_eqn = 3

    # Number of auxiliary variables in the aux array (initialized in setaux)
    clawdata.num_aux = 3
    
    # Index of aux array corresponding to capacity function, if there is one:
    clawdata.capa_index = 2
    
    
    # -------------
    # Initial time:
    # -------------

    clawdata.t0 = 0.0
    

    # Restart from checkpoint file of a previous run?
    # Note: If restarting, you must also change the Makefile to set:
    #    RESTART = True
    # If restarting, t0 above should be from original run, and the
    # restart_file 'fort.chkNNNNN' specified below should be in 
    # the OUTDIR indicated in Makefile.

    clawdata.restart = False                # True to restart from prior results
    clawdata.restart_file = 'fort.chk07323'  # File to use for restart data
    
    
    # -------------
    # Output times:
    #--------------

    # Specify at what times the results should be written to fort.q files.
    # Note that the time integration stops after the final output time.
 
    clawdata.output_style = 1
 
    if clawdata.output_style==1:
        # Output ntimes frames at equally spaced times up to tfinal:
        # Can specify num_output_times = 0 for no output
        clawdata.num_output_times = 20
        clawdata.tfinal = 20*3600
        clawdata.output_t0 = True  # output at initial (or restart) time?
        
    elif clawdata.output_style == 2:
        # Specify a list or numpy array of output times:
        # Include t0 if you want output at the initial time.
        #clawdata.output_times =  list(np.linspace(15*3600,32*3600,18))
        clawdata.output_times =  list(np.linspace(31*3600,60*3600,30))
        #clawdata.output_times =  list(np.linspace(31*3600,32*3600,13))
        #clawdata.output_times =  [14.5*3600] + list(np.linspace(15*3600,22*3600,85))
 
    elif clawdata.output_style == 3:
        # Output every step_interval timesteps over total_steps timesteps:
        clawdata.output_step_interval = 1
        clawdata.total_steps = 7323+10
        clawdata.output_t0 = True  # output at initial (or restart) time?
        

    clawdata.output_format = 'ascii'      # 'ascii', 'binary', 'netcdf'

    clawdata.output_q_components = 'all'   # could be list such as [True,True]
    clawdata.output_aux_components = 'none'  # could be list
    clawdata.output_aux_onlyonce = True    # output aux arrays only at t0
    

    # ---------------------------------------------------
    # Verbosity of messages to screen during integration:  
    # ---------------------------------------------------

    # The current t, dt, and cfl will be printed every time step
    # at AMR levels <= verbosity.  Set verbosity = 0 for no printing.
    #   (E.g. verbosity == 2 means print only on levels 1 and 2.)
    clawdata.verbosity = 0
    
    

    # --------------
    # Time stepping:
    # --------------

    # if dt_variable==True:  variable time steps used based on cfl_desired,
    # if dt_variable==Falseixed time steps dt = dt_initial always used.
    clawdata.dt_variable = True
    
    # Initial time step for variable dt.  
    # (If dt_variable==0 then dt=dt_initial for all steps)
    clawdata.dt_initial = 1
    
    # Max time step to be allowed if variable dt used:
    clawdata.dt_max = 1e+99
    
    # Desired Courant number if variable dt used 
    clawdata.cfl_desired = 0.5
    # max Courant number to allow without retaking step with a smaller dt:
    clawdata.cfl_max = 1.0
    
    # Maximum number of time steps to allow between output times:
    clawdata.steps_max = 50000


    # ------------------
    # Method to be used:
    # ------------------

    # Order of accuracy:  1 => Godunov,  2 => Lax-Wendroff plus limiters
    clawdata.order = 2
    
    # Use dimensional splitting? (not yet available for AMR)
    clawdata.dimensional_split = 'unsplit'
    
    # For unsplit method, transverse_waves can be 
    #  0 or 'none'      ==> donor cell (only normal solver used)
    #  1 or 'increment' ==> corner transport of waves
    #  2 or 'all'       ==> corner transport of 2nd order corrections too
    clawdata.transverse_waves = 2
    
    
    # Number of waves in the Riemann solution:
    clawdata.num_waves = 3
    
    # List of limiters to use for each wave family:  
    # Required:  len(limiter) == num_waves
    # Some options:
    #   0 or 'none'     ==> no limiter (Lax-Wendroff)
    #   1 or 'minmod'   ==> minmod
    #   2 or 'superbee' ==> superbee
    #   3 or 'vanleer'  ==> van Leer
    #   4 or 'mc'       ==> MC limiter
    clawdata.limiter = ['vanleer', 'vanleer', 'vanleer']
    
    clawdata.use_fwaves = True    # True ==> use f-wave version of algorithms
    
    # Source terms splitting:
    #   src_split == 0 or 'none'    ==> no source term (src routine never called)
    #   src_split == 1 or 'godunov' ==> Godunov (1st order) splitting used, 
    #   src_split == 2 or 'strang'  ==> Strang (2nd order) splitting used,  not recommended.
    clawdata.source_split = 1
    
    
    # --------------------
    # Boundary conditions:
    # --------------------

    # Number of ghost cells (usually 2)
    clawdata.num_ghost = 2
    
    # Choice of BCs at xlower and xupper:
    #   0 or 'user'     => user specified (must modify bcNamr.f to use this option)
    #   1 or 'extrap'   => extrapolation (non-reflecting outflow)
    #   2 or 'periodic' => periodic (must specify this at both boundaries)
    #   3 or 'wall'     => solid wall for systems where q(2) is normal velocity
    
    clawdata.bc_lower[0] = 'wall'   # at xlower
    clawdata.bc_upper[0] = 'wall'   # at xupper
    clawdata.bc_lower[1] = 'wall'   # at ylower
    clawdata.bc_upper[1] = 'wall'   # at yupper

    if 0:
        clawdata.bc_lower[0] = 'extrap'   # at xlower
        clawdata.bc_upper[0] = 'extrap'   # at xupper
        clawdata.bc_lower[1] = 'extrap'   # at ylower
        clawdata.bc_upper[1] = 'extrap'   # at yupper
                  
       
    # ---------------
    # Gauges:
    # ---------------

                  
    # --------------
    # Checkpointing:
    # --------------

    # Specify when checkpoint files should be created that can be
    # used to restart a computation.

    clawdata.checkpt_style = 1

    if clawdata.checkpt_style == 0:
      # Do not checkpoint at all
      pass

    elif clawdata.checkpt_style == 1:
      # Checkpoint only at tfinal.
      pass

    elif clawdata.checkpt_style == 2:
      # Specify a list of checkpoint times.  
      clawdata.checkpt_times = [30*3600.]  #list(3600.*np.array([12,24,36]))

    elif clawdata.checkpt_style == 3:
      # Checkpoint every checkpt_interval timesteps (on Level 1)
      # and at the final time.
      clawdata.checkpt_interval = 5

    

    # ---------------
    # AMR parameters:   (written to amr.data)
    # ---------------
    amrdata = rundata.amrdata

    # max number of refinement levels:
    amrdata.amr_levels_max = 1

    # List of refinement ratios at each level (length at least amr_level_max-1)
    amrdata.refinement_ratios_x = [6,4,3,30]
    amrdata.refinement_ratios_y = [6,4,3,30]
    amrdata.refinement_ratios_t = [6,4,3,30]  # not used


    # Specify type of each aux variable in amrdata.auxtype.
    # This must be a list of length num_aux, each element of which is one of:
    #   'center',  'capacity', 'xleft', or 'yleft'  (see documentation).
    amrdata.aux_type = ['center', 'capacity', 'yleft']


    # Flag for refinement based on Richardson error estimater:
    amrdata.flag_richardson = False    # use Richardson?
    amrdata.flag_richardson_tol = 1.0  # Richardson tolerance
    
    # Flag for refinement using routine flag2refine:
    amrdata.flag2refine = True      # use this?
    amrdata.flag2refine_tol = 0.5  # tolerance used in this routine
    # Note: in geoclaw the refinement tolerance is set as wave_tolerance below 
    # and flag2refine_tol is unused!

    # steps to take on each level L between regriddings of level L+1:
    amrdata.regrid_interval = 3       

    # width of buffer zone around flagged points:
    # (typically the same as regrid_interval so waves don't escape):
    amrdata.regrid_buffer_width  = 2

    # clustering alg. cutoff for (# flagged pts) / (total # of cells refined)
    # (closer to 1.0 => more small grids may be needed to cover flagged cells)
    amrdata.clustering_cutoff = 0.7

    # print info about each regridding up to this level:
    amrdata.verbosity_regrid = 0      


    # ---------------
    # Regions:
    # ---------------
    regions = rundata.regiondata.regions 
    # to specify regions of refinement append lines of the form
    #  [minlevel,maxlevel,t1,t2,x1,x2,y1,y2]

    regions.append([1, 2, 0., 1e9, 0, 360, -90, 90])                #whole world
    
    #  ----- For developers ----- 
    # Toggle debugging print statements:
    amrdata.dprint = False      # print domain flags
    amrdata.eprint = False      # print err est flags
    amrdata.edebug = False      # even more err est flags
    amrdata.gprint = False      # grid bisection/clustering
    amrdata.nprint = False      # proper nesting output
    amrdata.pprint = False      # proj. of tagged points
    amrdata.rprint = False      # print regridding summary
    amrdata.sprint = False      # space/memory output
    amrdata.tprint = False      # time step reporting each level
    amrdata.uprint = False      # update/upbnd reporting
    
    return rundata

    # end of function setrun
    # ----------------------


#-------------------
def setgeo(rundata):
#-------------------
    """
    Set GeoClaw specific runtime parameters.
    """

    try:
        geo_data = rundata.geo_data
    except:
        print "*** Error, this rundata has no geo_data attribute"
        raise AttributeError("Missing geo_data attribute")

    # == Physics ==
    geo_data.gravity = 9.81
    geo_data.coordinate_system =  2
    geo_data.earth_radius = 6367500.0

    # == Forcing Options
    geo_data.coriolis_forcing = False

    # == Algorithm and Initial Conditions ==
    #tide_stage = 0.
    #geo_data.sea_level = (tide_stage - 77.)/100.    #  m relative to MHW
    geo_data.sea_level = 0.

    geo_data.dry_tolerance = 0.001
    geo_data.friction_forcing = True
    geo_data.manning_coefficient = 0.025
    geo_data.friction_depth = 100.0

    # Refinement settings
    refinement_data = rundata.refinement_data
    refinement_data.variable_dt_refinement_ratios = True
    refinement_data.wave_tolerance = 0.02
    refinement_data.deep_depth = 100.0
    refinement_data.max_level_deep = 4

    # == settopo.data values ==
    topofiles = rundata.topo_data.topofiles
    # for topography, append lines of the form
    #    [topotype, minlevel, maxlevel, t1, t2, fname]

    try:
        #CCdir = os.environ['CCdir']
        CCdir = '.'
    except: 
        raise Exception("You must set the environment variable CCdir")

    # BELOW FOR CHILE RUN
    #new topo file in CCptha2:
    topofiles.append([3, 1, 1, 0., 1.e10, 'small_topo.tt3'])

                              #CCdir + '/topo/etopo4min120E65W65S65N.asc'])


    # == setdtopo.data values ==
    dtopofiles = rundata.dtopo_data.dtopofiles
    # for moving topography, append lines of the form :  
    #   [topotype, minlevel,maxlevel,fname]
    #CCdir = os.environ['CCdir'] 
    #dtopofile = CCdir + '/dtopo/SChSZ/SChSZe14v2.tt3'
    #dtopotype  = 3
    #dtopofiles.append([3,1,4,dtopofile])
    #rundata.dtopo_data.dt_max_dtopo = 0.45 # dt while topo is moving


    # == setqinit.data values ==
    rundata.qinit_data.qinit_type =  4
    rundata.qinit_data.qinitfiles = []
    qinitfiles = rundata.qinit_data.qinitfiles 
    # for qinit perturbations, append lines of the form: (<= 1 allowed for now!)
    #   [minlev, maxlev, fname]
    rundata.qinit_data.qinitfiles.append([1, 2, 'hump.xyz'])


    # == fixedgrids.data values ==
    rundata.fixed_grid_data.fixedgrids = []
    fixedgrids = rundata.fixed_grid_data.fixedgrids
    # for fixed grids append lines of the form
    # [t1,t2,noutput,x1,x2,y1,y2,xpoints,ypoints,\
    #  ioutarrivaltimes,ioutsurfacemax]

    # == fgmax.data values ==
    fgmax_files = rundata.fgmax_data.fgmax_files
    # for fixed grids append to this list names of any fgmax input files
    #fgmax_files.append('fgmax_grid1.txt')  
    #fgmax_files.append('fgmax_grid2.txt')  
    #rundata.fgmax_data.num_fgmax_val = 5



    return rundata
    # end of function setgeo
    # ----------------------


from clawpack.geoclaw import fgmax_tools as F
from numpy import floor, ceil

tstart_max = max(100., t_harbor)     # after dtopo motion

def make_fgmax_grid1():

    FG = F.fgmax_grid_parameters()
    FG.fname = 'fgmax_grid1.txt'
    FG.point_style = 2       # will specify a 2d grid of points

    # rough outline of desired fgmax grid (will be adjusted below):
    cellsize = 9.259260e-05   # cellsize  from cc_1_3 topo file
    xlower = 235.77 + cellsize
    xupper = 235.85
    ylower = 41.735 + cellsize
    yupper = 41.785

    FG.tstart_max =  tstart_max  # when to start monitoring max values
    FG.tend_max = 1.e10       # when to stop monitoring max values
    FG.dt_check = 30.       # target time increment between updating max values
    FG.min_level_check = 6   # which levels to monitor max on
    FG.arrival_tol = 1.e-2        # tolerance for flagging arrival

    # what's below is set up for CC 1/3" grid:

    # Make sure fgmax points align exactly with 1/3" topo grid points:
    xlow_topo = 2.357655e+02  # xll from cc_1_3 topo file
    ylow_topo = 4.171643e+01  # yll  from cc_1_3 topo file
    cellsize = 9.259260e-05   # cellsize  from cc_1_3 topo file
    #dx = 1./(3*3600.)    # 1/3 arcssecond grid

    FG.dx = cellsize 

    ilower = int(floor((xlower - xlow_topo)/cellsize + 1e-6))
    iupper = int(ceil((xupper - xlow_topo)/cellsize - 1e-6))
    mx = iupper - ilower + 1
    FG.x1 = xlow_topo + ilower*cellsize
    FG.x2 = xlow_topo + iupper*cellsize
    jlower = int(floor((ylower - ylow_topo)/cellsize + 1e-6))
    jupper = int(ceil((yupper - ylow_topo)/cellsize - 1e-6))
    my = jupper - jlower + 1
    FG.y1 = ylow_topo + jlower*cellsize
    FG.y2 = ylow_topo + jupper*cellsize
    
    print "Adjusted corners of grid to match topo file: "
    print "  x1 = %15.10f, y1 = %15.10f" % (FG.x1,FG.y1)
    print "  x2 = %15.10f, y2 = %15.10f" % (FG.x2,FG.y2)
    print "  cellsize = %g" % cellsize

    F.make_fgmax(FG)


if __name__ == '__main__':
    # Set up run-time parameters and write all data files.
    import sys
    rundata = setrun(*sys.argv[1:])
    rundata.write()

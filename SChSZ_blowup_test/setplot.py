
""" 
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.
    
""" 


from clawpack.geoclaw import topotools
import pylab
import glob
from numpy import loadtxt, linspace



# --------------------------
def setplot(plotdata):
# --------------------------
    
    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of pyclaw.plotters.data.ClawPlotData.
    Output: a modified version of plotdata.
    
    """ 


    from clawpack.visclaw import colormaps
    from clawpack.geoclaw import geoplot

    plotdata.clearfigures()  # clear any old figures,axes,items dat
    plotdata.format = "ascii"

    try:
        tsudata = open(plotdata.outdir+'/geoclaw.data').readlines()
        for line in tsudata:
            if 'sea_level' in line:
                sea_level = float(line.split()[0])
                print "sea_level = ",sea_level
    except:
        print "Could not read sea_level, setting to 0."
        sea_level = 0.


    clim_ocean = 0.5
    clim_CC = 0.5

    cmax_ocean = clim_ocean + sea_level
    cmin_ocean = -clim_ocean + sea_level
    cmax_CC = clim_CC + sea_level
    cmin_CC = -clim_CC + sea_level


    # To plot gauge locations on pcolor or contour plot, use this as
    # an afteraxis function:

    def addgauges(current_data):
        from clawpack.visclaw import gaugetools
        gaugetools.plot_gauge_locations(current_data.plotdata, \
             gaugenos='all', format_string='ko', add_labels=True)

    def timeformat(t):
        from numpy import mod
        hours = int(t/3600.)
        tmin = mod(t,3600.)
        min = int(tmin/60.)
        sec = int(mod(tmin,60.))
        timestr = '%s:%s:%s' % (hours,str(min).zfill(2),str(sec).zfill(2))
        return timestr
        
    def title_hours(current_data):
        from pylab import title
        t = current_data.t
        timestr = timeformat(t)
        title('%s after earthquake' % timestr)

    def plotcc(current_data):
        from pylab import plot,text
        plot([235.8162], [41.745616],'wo')
        text(235.8,41.9,'Cr.City',color='w',fontsize=10)
    

    #-----------------------------------------
    # Figure for big area
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='Pacific', figno=0)
    #plotfigure.kwargs = {'figsize': (16,4)}
    #plotfigure.show = False

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'Pacific'
    plotaxes.scaled = True

    def aa(current_data):
        from pylab import ticklabel_format, xticks, gca, cos, pi
        #plotcc(current_data)
        title_hours(current_data)
        ticklabel_format(format='plain',useOffset=False)
        xticks(rotation=20)
        a = gca()
        a.set_aspect(1./cos(41.75*pi/180.))
    plotaxes.afteraxes = aa

    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
    plotitem.plot_var = geoplot.surface_or_depth
    my_cmap = colormaps.make_colormap({-1.0: [0.0,0.0,1.0], \
                                     -0.5: [0.5,0.5,1.0], \
                                      0.0: [1.0,1.0,1.0], \
                                      0.5: [1.0,0.5,0.5], \
                                      1.0: [1.0,0.0,0.0]})
    plotitem.imshow_cmap = my_cmap
    plotitem.imshow_cmin = -0.5
    plotitem.imshow_cmax = 0.5
    plotitem.add_colorbar = True
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.amr_patchedges_show = [0]

    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
    plotitem.plot_var = geoplot.land
    plotitem.imshow_cmap = geoplot.land_colors
    plotitem.imshow_cmin = 0.0
    plotitem.imshow_cmax = 100.0
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.amr_patchedges_show = [0]
    #plotaxes.xlimits = [228,238] 
    #plotaxes.ylimits = [34,50]
    #plotaxes.afteraxes = addgauges

    # Add contour lines of bathymetry:
    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    plotitem.show = False
    plotitem.plot_var = geoplot.topo
    from numpy import arange, linspace
    plotitem.contour_levels = linspace(-6000,0,7)
    plotitem.amr_contour_colors = ['g']  # color on each level
    plotitem.kwargs = {'linestyles':'solid'}
    plotitem.amr_contour_show = [0,0,1,0]  # show contours only on finest level
    plotitem.celledges_show = 0
    plotitem.patchedges_show = 0

    # Add contour lines of topography:
    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    plotitem.show = False
    plotitem.plot_var = geoplot.topo
    from numpy import arange, linspace
    plotitem.contour_levels = arange(0., 11., 1.)
    plotitem.amr_contour_colors = ['g']  # color on each level
    plotitem.kwargs = {'linestyles':'solid'}
    plotitem.amr_contour_show = [0,0,0,1]  # show contours only on finest level
    plotitem.celledges_show = 0
    plotitem.patchedges_show = 0


    # =====

    plotfigure = plotdata.new_plotfigure(name='Tahiti', figno=20)
    #plotfigure.kwargs = {'figsize': (16,4)}
    #plotfigure.show = False

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'Pacific'
    plotaxes.scaled = True

    def aa(current_data):
        from pylab import ticklabel_format, xticks, gca, cos, pi
        #plotcc(current_data)
        title_hours(current_data)
        ticklabel_format(format='plain',useOffset=False)
        xticks(rotation=20)
        a = gca()
        a.set_aspect(1./cos(41.75*pi/180.))
    plotaxes.afteraxes = aa

    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
    plotitem.plot_var = geoplot.surface_or_depth
    my_cmap = colormaps.make_colormap({-1.0: [0.0,0.0,1.0], \
                                     -0.5: [0.5,0.5,1.0], \
                                      0.0: [1.0,1.0,1.0], \
                                      0.5: [1.0,0.5,0.5], \
                                      1.0: [1.0,0.0,0.0]})
    plotitem.imshow_cmap = my_cmap
    plotitem.imshow_cmin = -0.5
    plotitem.imshow_cmax = 0.5
    plotitem.add_colorbar = True
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.amr_patchedges_show = [0]

    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
    plotitem.plot_var = geoplot.land
    plotitem.imshow_cmap = geoplot.land_colors
    plotitem.imshow_cmin = 0.0
    plotitem.imshow_cmax = 100.0
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.amr_patchedges_show = [0]
    plotaxes.xlimits = [215.5,217.5] 
    plotaxes.ylimits = [-18,-15.5]
    #plotaxes.afteraxes = addgauges

    # =====
    # Transects:
    latitude = -16.65
    longitude = 216.42

    plotfigure = plotdata.new_plotfigure(name='x-transect', figno=41)
    #plotfigure.show = False
    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(211)'
    plotaxes.title = 'Transect in x at latitude y = %s' % latitude

    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')

    def slice_x(current_data):
        from pylab import find,ravel
        x = current_data.x
        y = current_data.y
        dy = current_data.dy
        q = current_data.q
        ij = find((y <= latitude+dy/2.) & (y > latitude-dy/2.))
        x_slice = ravel(x)[ij]
        eta_slice = ravel(q[3,:,:])[ij]
        print "+++ min eta = ",eta_slice.min()
        return x_slice, eta_slice

    plotitem.map_2d_to_1d = slice_x

    def slice_xp(current_data):
        from pylab import find,ravel
        x = current_data.x
        y = current_data.y
        dy = current_data.dy
        q = current_data.q
        ij = find((y <= latitude+dy+dy/2.) & (y > latitude+dy-dy/2.))
        x_slice = ravel(x)[ij]
        eta_slice = ravel(q[3,:,:])[ij]
        return x_slice, eta_slice

    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')
    plotitem.map_2d_to_1d = slice_xp
    plotitem.color = 'r'


    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(212)'
    plotaxes.title = ''



    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')
    plotitem.map_2d_to_1d = slice_x
    plotitem.color = 'b'

    def B_x(current_data):
        from pylab import find,ravel
        x = current_data.x
        y = current_data.y
        dy = current_data.dy
        q = current_data.q
        ij = find((y <= latitude+dy/2.) & (y > latitude-dy/2.))
        x_slice = ravel(x)[ij]
        eta_slice = ravel(q[3,:,:])[ij]
        h_slice = ravel(q[0,:,:])[ij]
        B_slice = eta_slice - h_slice
        return x_slice, B_slice

    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')
    plotitem.map_2d_to_1d = B_x
    plotitem.color = 'g'

    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')
    plotitem.map_2d_to_1d = slice_xp
    plotitem.color = 'r'

    def B_x(current_data):
        from pylab import find,ravel
        x = current_data.x
        y = current_data.y
        dy = current_data.dy
        q = current_data.q
        ij = find((y <= latitude+dy+dy/2.) & (y > latitude+dy-dy/2.))
        x_slice = ravel(x)[ij]
        eta_slice = ravel(q[3,:,:])[ij]
        h_slice = ravel(q[0,:,:])[ij]
        B_slice = eta_slice - h_slice
        return x_slice, B_slice

    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')
    plotitem.map_2d_to_1d = B_x
    plotitem.color = 'r'
        
    plotfigure = plotdata.new_plotfigure(name='y-transect', figno=22)
    #plotfigure.show = False
    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(211)'
    plotaxes.title = 'Transect in y at longitude x = %s' % longitude

    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')

    def slice_y(current_data):
        from pylab import find,ravel
        x = current_data.x
        y = current_data.y
        dx = current_data.dx
        q = current_data.q
        ij = find((x <= longitude+dx/2.) & (x > longitude-dx/2.))
        y_slice = ravel(y)[ij]
        eta_slice = ravel(q[3,:,:])[ij]
        return y_slice, eta_slice

    plotitem.map_2d_to_1d = slice_y

    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(212)'
    plotaxes.title = ''
    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')
    plotitem.map_2d_to_1d = slice_y
    plotitem.color = 'b'

    def B_y(current_data):
        from pylab import find,ravel
        x = current_data.x
        y = current_data.y
        dy = current_data.dy
        q = current_data.q
        dx = current_data.dx
        q = current_data.q
        ij = find((x <= longitude+dx/2.) & (x > longitude-dx/2.))
        y_slice = ravel(y)[ij]
        eta_slice = ravel(q[3,:,:])[ij]
        h_slice = ravel(q[0,:,:])[ij]
        B_slice = eta_slice - h_slice
        return y_slice, B_slice

    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')
    plotitem.map_2d_to_1d = B_y
    plotitem.color = 'g'
        
        
    plotfigure = plotdata.new_plotfigure(name='x-velocity', figno=23)
    #plotfigure.show = False
    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(211)'
    plotaxes.title = 'u in x at latitude y = %s' % latitude

    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')

    def slice_x(current_data):
        from pylab import find,ravel,where
        x = current_data.x
        y = current_data.y
        dy = current_data.dy
        q = current_data.q
        ij = find((y <= latitude+dy/2.) & (y > latitude-dy/2.))
        x_slice = ravel(x)[ij]
        eta_slice = ravel(q[3,:,:])[ij]
        h_slice = ravel(q[0,:,:])[ij]
        hu_slice = ravel(q[1,:,:])[ij]
        u_slice = where(h_slice>0, hu_slice/h_slice, 0.)
        return x_slice, u_slice

    plotitem.map_2d_to_1d = slice_x

    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(212)'
    plotaxes.title = 'hu'


    def slice_x(current_data):
        from pylab import find,ravel
        x = current_data.x
        y = current_data.y
        dy = current_data.dy
        q = current_data.q
        ij = find((y <= latitude+dy/2.) & (y > latitude-dy/2.))
        x_slice = ravel(x)[ij]
        eta_slice = ravel(q[3,:,:])[ij]
        h_slice = ravel(q[0,:,:])[ij]
        hu_slice = ravel(q[1,:,:])[ij]
        return x_slice, hu_slice

    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')
    plotitem.map_2d_to_1d = slice_x
    plotitem.color = 'b'
        

    plotfigure = plotdata.new_plotfigure(name='v-transect', figno=24)
    #plotfigure.show = False
    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(211)'
    plotaxes.title = 'v in y at longitude x = %s' % longitude

    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')

    def slice_y(current_data):
        from pylab import find,ravel,where
        x = current_data.x
        y = current_data.y
        dx = current_data.dx
        q = current_data.q
        ij = find((x <= longitude+dx/2.) & (x > longitude-dx/2.))
        y_slice = ravel(y)[ij]
        eta_slice = ravel(q[3,:,:])[ij]
        h_slice = ravel(q[0,:,:])[ij]
        hv_slice = ravel(q[2,:,:])[ij]
        v_slice = where(h_slice>0, hv_slice/h_slice, 0.)
        return y_slice, v_slice

    plotitem.map_2d_to_1d = slice_y

    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(212)'
    plotaxes.title = 'hv'


    def slice_y(current_data):
        from pylab import find,ravel,where
        x = current_data.x
        y = current_data.y
        dx = current_data.dx
        q = current_data.q
        ij = find((x <= longitude+dx/2.) & (x > longitude-dx/2.))
        y_slice = ravel(y)[ij]
        eta_slice = ravel(q[3,:,:])[ij]
        h_slice = ravel(q[0,:,:])[ij]
        hv_slice = ravel(q[2,:,:])[ij]
        v_slice = where(h_slice>0, hv_slice/h_slice, 0.)
        return y_slice, hv_slice

        
    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')
    plotitem.map_2d_to_1d = slice_y
    plotitem.color = 'b'


    #-----------------------------------------
    # Figure for zoom
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='California', figno=10)
    plotfigure.show = False

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'California'
    plotaxes.scaled = True
    plotaxes.afteraxes = aa

    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
    plotitem.plot_var = geoplot.surface_or_depth
    plotitem.imshow_cmap = my_cmap
    plotitem.imshow_cmin = cmin_ocean
    plotitem.imshow_cmax = cmax_ocean
    plotitem.add_colorbar = True
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.amr_patchedges_show = [1]

    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
    plotitem.plot_var = geoplot.land
    plotitem.imshow_cmap = geoplot.land_colors
    plotitem.imshow_cmin = 0.0
    plotitem.imshow_cmax = 100.0
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.amr_patchedges_show = [1]
    plotaxes.xlimits = [232,237] 
    plotaxes.ylimits = [39,44]
    plotaxes.afteraxes = aa

    
 
    #-----------------------------------------
    # Figure for zoom2
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='Crescent City', figno=11)
    plotfigure.show = False

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'Crescent City'
    plotaxes.scaled = True

    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
    plotitem.plot_var = geoplot.surface_or_depth
    plotitem.imshow_cmap = my_cmap
    plotitem.imshow_cmin = cmin_CC
    plotitem.imshow_cmax = cmax_CC
    plotitem.add_colorbar = True
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.amr_patchedges_show = [1]

    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
    plotitem.plot_var = geoplot.land
    plotitem.imshow_cmap = geoplot.land_colors
    plotitem.imshow_cmin = 0.0
    plotitem.imshow_cmax = 100.0
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.amr_patchedges_show = [1]
    plotaxes.xlimits = [235.76,235.84] 
    plotaxes.ylimits = [41.72,41.78]
    def aa(current_data):
        from pylab import ticklabel_format, xticks, gca, cos, pi
        addgauges(current_data)
        title_hours(current_data)
        ticklabel_format(format='plain',useOffset=False)
        xticks(rotation=20)
        a = gca()
        a.set_aspect(1./cos(41.75*pi/180.))
    plotaxes.afteraxes = aa

    # add contour lines of bathy if desired:
    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    #plotitem.show = False
    plotitem.plot_var = geoplot.topo
    plotitem.contour_levels = [0.]
    plotitem.amr_contour_colors = ['k']  # color on each level
    plotitem.kwargs = {'linestyles':'solid','linewidths':2}
    plotitem.amr_contour_show = [0,0,0,0,0,1]
    plotitem.celledges_show = 0
    plotitem.patchedges_show = 0

 

    #-----------------------------------------
    # Figures for gauges
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='gauge plot', figno=300, \
                    type='each_gauge')
    #plotfigure.clf_each_gauge = False

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'Surface'

    # Plot surface as blue curve:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 3
    plotitem.plotstyle = 'b-'


    def fix_gauge(current_data):
        from pylab import plot, legend, xticks, floor, yticks,xlabel,savefig
        t = current_data.t
        gaugeno = current_data.gaugeno
        n = int(floor(t.max()/1800.) + 2)
        xticks([1800*i for i in range(n)],[str(i/2.) for i in range(n)],\
          fontsize=15)
        yticks(fontsize=15)
        xlabel("Hours")

    #plotaxes.afteraxes = fix_gauge




    #-----------------------------------------
    
    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via pyclaw.plotters.frametools.printframes:

    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos = 'all'  # range(0,50,2)
    plotdata.print_fignos = 'all'            # list of figures to print
    plotdata.print_gaugenos = 'all'          # list of gauges to print
    plotdata.html = True                     # create html files of plots?
    plotdata.html_homelink = '../README.html'   # pointer for top of index
    plotdata.latex = True                    # create latex file of plots?
    plotdata.latex_figsperline = 2           # layout of plots
    plotdata.latex_framesperline = 1         # layout of plots
    plotdata.latex_makepdf = False           # also run pdflatex?

    return plotdata

    

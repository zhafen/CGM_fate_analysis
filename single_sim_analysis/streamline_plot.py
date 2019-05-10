#!/usr/bin/env python
# coding: utf-8

# In[1]:


import copy
import numpy as np
import sys


# In[2]:


import analysis_config


# In[3]:


import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.cm as cm


# In[4]:


import linefinder.analyze_data.worldlines as a_worldlines
import linefinder.analyze_data.plot_worldlines as p_worldlines
import linefinder.utils.presentation_constants as p_constants


# In[5]:


import galaxy_dive.analyze_data.ahf as analyze_ahf
import galaxy_dive.plot_data.ahf as plot_ahf
import galaxy_dive.analyze_data.particle_data as particle_data
import galaxy_dive.plot_data.generic_plotter as generic_plotter
import galaxy_dive.plot_data.plotting as plotting
import galaxy_dive.utils.data_operations as data_operations


# In[6]:


import linefinder.utils.file_management as file_management
import linefinder.config as config


# # Load Data

# In[7]:


default_sim = 'm12i'
default_snum = 465


# In[8]:


if len( sys.argv ) == 2:
    sim_name = sys.argv[1]
else:
    sim_name = default_sim


# In[9]:


# Automatically set or retrieve args
if len( sys.argv ) == 3:
    try:
        snum = int( sys.argv[2] )
    except ValueError:
        snum = default_snum
else:
    snum = default_snum


# In[10]:


galdef = ''


# In[11]:


file_manager = file_management.FileManager( 'CGM_fate' )


# In[12]:


file_manager.get_linefinder_dir( sim_name, )


# In[13]:


defaults = file_manager.get_linefinder_analysis_defaults(
    '_CGM_snum{}'.format( snum ),
    sim_name = sim_name,
    galdef = galdef
)


# In[14]:


ind = defaults['ahf_index'] - snum


# In[15]:


w = a_worldlines.Worldlines( **defaults )


# In[16]:


w.retrieve_halo_data()


# In[17]:


m_plot_label  = r'$M_{\rm h} = 10^{' + '{:.02g}'.format( np.log10( w.m_vir[snum] ) )
m_plot_label += '} M_\odot$'
plot_label = m_plot_label + ', z={:.02}'.format( w.redshift[snum] )
print( plot_label )


# In[18]:


classification_list = copy.copy( p_constants.CLASSIFICATIONS_CGM_FATE )


# In[19]:


w_plotter = p_worldlines.WorldlinesPlotter( w, label=plot_label )


# In[20]:


g_data = particle_data.ParticleData(
    sdir = file_manager.get_sim_dir( sim_name ),
    snum = snum,
    ptype = config.PTYPE_GAS,
    halo_data_dir = file_manager.get_halo_dir( sim_name ),
    main_halo_id = config.MAIN_MT_HALO_ID[sim_name],
)


# In[21]:


g_plotter = generic_plotter.GenericPlotter(
    g_data,
    label=plot_label,
)


# In[22]:


s_data = particle_data.ParticleData(
    sdir = file_manager.get_sim_dir( sim_name ),
    snum = snum,
    ptype = config.PTYPE_STAR,
    halo_data_dir = file_manager.get_halo_dir( sim_name ),
    main_halo_id = config.MAIN_MT_HALO_ID[sim_name],    
)


# In[23]:


s_plotter = generic_plotter.GenericPlotter( s_data )


# ### Create a circle to plot

# In[24]:


w.halo_data.data_dir


# In[25]:


r_gal = w.r_gal[ind]


# In[26]:


circle = []
for phi in np.linspace( 0., 2.*np.pi, 256 ):
    
    circle.append(
        [ r_gal*np.cos(phi), r_gal*np.sin(phi), 0. ]
    )
    
circle = np.array( circle )

rotated_circle = data_operations.align_axes( circle, s_data.total_ang_momentum, )


# # Illustrative Plot

# In[27]:


r_vir = w.r_vir.values[ind]


# In[28]:


t_show_min = {
    465 : 0.5,
    172 : 0.25,
    214 : 0.25,
}
t_show_max = {
    465 : 1.0,
    172 : 0.5,
    214 : 0.5,
}


# In[29]:


data_args = {
    465 : { 'smooth_data' : True, 'smoothing_window_length' : 21 },
    172 : { 'smooth_data' : True, 'smoothing_window_length' : 21 },
    214 : { 'smooth_data' : True, 'smoothing_window_length' : 21 },
}


# In[30]:


fig = plt.figure( figsize=(20,10), facecolor='white' )
main_ax = plt.gca()

gs = matplotlib.gridspec.GridSpec( 1, 2 )
axs = [ plt.subplot( gs[0,0] ), plt.subplot( gs[0,1] ) ]

gs.update( wspace=0.15 )

plotted_range = [ -1.2*r_vir, 1.2*r_vir ]

x_axis = 'Rx'
y_axis = 'Ry'

for j, sample_selected_interval in enumerate( [ False, True ] ):
    
    ax = axs[j]
    
    labels = []
    color_objects = []
    for i, classification in enumerate( classification_list ):

        w_plotter.data_object.data_masker.mask_data( 'PType', data_value=config.PTYPE_GAS )
        
        # Don't plot CGM EP for m10s at low-z when applying the time cut, because there are just
        # too few particles
        if (
            classification == 'is_outside_any_gal_EP' and
            ax.is_last_col() and
            config.MASS_BINS[sim_name] == 'm10' and
            snum == 465
        ):
            continue
            
        # Don't plot halo transfer, as it's an unimportant component
        if classification == 'is_CGM_halo_transfer':
            continue
            
        # Don't plot CGM still when we're doing fates close to their exit time
        if classification == 'is_CGM_still' and sample_selected_interval:
            continue
        
        print( classification )
        w_plotter.plot_streamlines(
            x_axis,
            y_axis,
            ax = ax,
            classification = classification,
            classification_ind = ind,
            start_ind = 'time_based',
            end_ind = ind,
            t_start = 1.0,
            sample_size = 500,
            sample_selected_interval = sample_selected_interval,
            selected_interval_type = 'time_until_not',
            selected_interval_classification = 'is_in_CGM_or_interface',
            t_show_min = t_show_min[snum],
            t_show_max = t_show_max[snum],
            x_data_kwargs = data_args[snum],
            y_data_kwargs = data_args[snum],
            linewidth = 2.5,
            plot_halos = False,
#             plot_xlabel = ( sim_name == 'm10y' ),
            plot_ylabel = ax.is_first_col(),
            x_label = '{} position (pkpc)'.format( x_axis[1] ),
            y_label = '{} position (pkpc)'.format( y_axis[1] ),
            x_range = plotted_range,
            y_range = plotted_range,
            fontsize = 24,
        )

        # Make virtual artists to allow a legend to appear
        color_object = matplotlib.patches.Rectangle(                         
            (0, 0),                                                          
            1,                                                               
            1,                                                               
            fc = p_constants.CLASSIFICATION_COLORS_B[classification],                                 
            ec = p_constants.CLASSIFICATION_COLORS_B[classification],                                 
            alpha = p_constants.CLASSIFICATION_ALPHA,                        
        )                                                                    
        color_objects.append( color_object )                                 
        labels.append( p_constants.CLASSIFICATION_LABELS[classification] )

    if j==0:
        plot_label = sim_name
    else:
        plot_label = r'$z = {:.2g}$'.format( g_data.redshift )
    vmaxes = {
        'm10' : 4e3,
        'm11' : 4e3,
        'm12' : 4e3,
    }
    vmax = vmaxes[config.MASS_BINS[sim_name]]
    g_plotter.histogram2d(
        x_axis,
        y_axis,
        cmap = cm.Greys,
        ax = ax,
        x_range = plotted_range,
        y_range = plotted_range,
        n_bins = 200,
#         vmin = 1,
        vmax = vmax,
        add_colorbar = False,
        add_x_label = False,
        add_y_label = False,
        label_fontsize = 24,
        plot_label = plot_label,
    )
    
    vmins = {
        172 : {
            'm10' : 10,
            'm11' : 10,
            'm12' : 20,
        },
        214 : {
            'm12' : 50,
        },
        465 : {
            'm10' : 20,
            'm11' : 20,
            'm12' : 5e2,
        },
    }
    vmin = vmins[snum][config.MASS_BINS[sim_name]]
    s_plotter.histogram2d(
        x_axis,
        y_axis,
        cmap = cm.viridis,
        ax = ax,
        x_range = plotted_range,
        y_range = plotted_range,
        n_bins = 200,
        vmin = vmin,
#         vmax = 4e3,
        plot_label = None,
        add_colorbar = False,
        add_x_label = False,
        add_y_label = False,
        label_fontsize = 20,
        zorder = 150,
        min_bin_value_displayed = vmin,
    )

#     s_plotter.scatter(
#         'Rx',
#         'Ry',
#         color = '#303030',
# #         color = '#e6e032',
#         marker = '*',
#         ax = ax,
#         x_range = plotted_range,
#         y_range = plotted_range,
#         n_subsample = 2000,
#         zorder = 150.,
#         plot_label = None,
#         add_x_label = False,
#         add_y_label = False,
#     )

    # Virial Radius circle
    cir = mpatches.Circle(
        [0, 0],
        radius = r_vir,
        linewidth = 3,
        color = 'w',
        linestyle = '--',
        fill = False,
        facecolor = 'w',
    )
    ax.add_patch( cir )
    ax.annotate(
        s = r'$R_{\rm vir}$',
        xy = ( .7475*r_vir, -.7475*r_vir ),
        xycoords = 'data',
        fontsize = 28,
        color = 'k',
        zorder = 150,
    )
    
    # Galaxy disk circle
    ax.plot(
        rotated_circle[:,0],
        rotated_circle[:,1],
        color = 'w',
        linewidth = 3,
        zorder = 200,
    )
    
    if ax.is_first_col():
        filter_label = 'No time selection'
    if ax.is_last_col():
        filter_label = 'Will leave CGM in {}-{} Gyr'.format(
            t_show_min[snum],
            t_show_max[snum],
        )
    ax.annotate(
        s = filter_label,
        xy = ( 1., 1. ),
        xycoords = 'axes fraction',
        va = 'bottom',
        ha = 'right',
        fontsize = 22,
    )

    if ax.is_first_col():
        leg = ax.legend(
            color_objects,
            labels,
            prop={'size': 20},
            ncol=1,
    #         loc=(0.65, 0.83),
            loc='upper left',
            fontsize=24,
            framealpha = 0.9,
        )
        leg.set_zorder( 200 )

    ax.set_aspect( 'equal' )

save_file = 'streamlines_{}.png'.format( defaults['tag'] )
plotting.save_fig(
    out_dir = file_manager.get_project_figure_dir(),
    save_file = save_file,
    fig = fig,
    resolution = 50,
)

fig


# ### Plot vs Time

# In[31]:


w.data_masker.clear_masks()


# In[32]:


y_max = np.nanpercentile( w.get_selected_data_over_time(
    data_key = 'R',
    snum = 465,
    classification = 'is_CGM_IGM_accretion',
)[:ind], 99. )*1.1


# In[40]:


used_classification_list = copy.copy( classification_list )
used_classification_list.remove( 'is_CGM_halo_transfer' )


# In[49]:


fig = plt.figure( figsize=(10,15), facecolor='white' )
main_ax = plt.gca()

gs = matplotlib.gridspec.GridSpec( 4, 1 )

y_maxes = []
for i, classification in enumerate( used_classification_list ):
    
    # Minor category, skip
    if classification == 'is_CGM_halo_transfer':
        continue
    
    ax = plt.subplot(gs[i,0])
    
#     if ax.is_first_row():
#         ax.annotate(
#             s = m_plot_label,
#             xy = (0,1.1),
#             xycoords = 'axes fraction',
#             va = 'bottom',
#             fontsize = 24,
#         )
    if classification == 'is_CGM_IP':
        ax.annotate(
            s = m_plot_label,
            xy = (0.05,0.95),
            xycoords = 'axes fraction',
            va = 'top',
            fontsize = 24,
        )
    
    gs.update(wspace=0.025, hspace=0.0001)
    
    w_plotter.data_object.data_masker.mask_data( 'PType', data_value=config.PTYPE_GAS )
    
    class_y_max = w_plotter.plot_streamlines_vs_time(
        y_key = 'R',
        classification = classification,
        classification_ind = ind,
        start_ind = 0,
        end_ind = 595,
        sample_size = 10,
        y_data_kwargs = { 'smooth_data' : True },
        ax = ax,
        x_range = [0.5, 13.8 ],
        y_range = [0., 550, ],
        horizontal_line_value = None,
        plot_CGM_region = True,
        CGM_region_alpha = 0.1,
        return_y_max = True,
#         line_features = {
#             'is_class': {
#                 'key': classification,
#                 'value': True,
#                 'data_kwargs': {},
#                 'line_attrs': {
#                     'linewidth': 3,
#                     'color': p_constants.CLASSIFICATION_COLORS_B[classification],
#                     },
#             },
#             'is_gas': {
#                 'key': classification,
#                 'value': False,
#                 'data_kwargs': {},
#                 'line_attrs': {
#                     'linewidth': 1,
#                     'color': p_constants.CLASSIFICATION_COLORS_B[classification],
#                     'linestyle': '--',
#                     },
#                 }
#         },
    )
    y_maxes.append( class_y_max )
        
    # Make virtual artists to allow a legend to appear
    labels = []
    color_objects = []
    color_object = matplotlib.patches.Rectangle(                         
        (0, 0),                                                          
        1,                                                               
        1,                                                               
        fc = p_constants.CLASSIFICATION_COLORS_B[classification],                                 
        ec = p_constants.CLASSIFICATION_COLORS_B[classification],                                 
        alpha = p_constants.CLASSIFICATION_ALPHA,                        
    )                                                                    
    color_objects.append( color_object )                                 
    labels.append( p_constants.CLASSIFICATION_LABELS[classification] )
    
    # Plot the outer edge of the galaxy under galdefv2
#     ax.plot( w.get_data( 'time')[:r_gal.size], w.r_gal, color='k', linewidth=3.5 )

    # Legend
    l = ax.legend(
        color_objects,
        labels,
        prop={'size': 18.5},
        ncol=1,
        bbox_to_anchor=(1, 0.95),
        loc='upper right',
        fontsize=20,
        framealpha=1.,
    )
    l.set_zorder( 200. )
        
# Changes to axes
for i, classification in enumerate( used_classification_list ):
    
    ax = plt.subplot(gs[i,0])
    ax.set_ylim( 0, 1.05*max( y_maxes ) )
    
    # Add redshift to the axis
    ax2 = plotting.add_redshift_to_axis(
        ax,
        hubble = w.ptracks.data_attrs['hubble'],
        omega_matter = w.ptracks.data_attrs['omega_matter'],
        tick_redshifts = np.array([ 0.1, 0.25, 0.4, 0.5, 0.75, 1, 1.5, 2, 3, 5, ])
    )
    
    # Change tick orientation
    ax.tick_params( direction='in' )
    ax2.tick_params( direction='in' )
    
    # Hide overlapping labels
    if not ax.is_first_row():
        ax2.xaxis.set_ticklabels([])
    if not ax.is_last_col():
        ax.xaxis.set_ticklabels([])
        
    # Avoid overlapping ticks
    ax.get_yticklabels()[0].set_verticalalignment( 'bottom' )
    ax.get_yticklabels()[1].set_verticalalignment( 'top' )
    
    # X label
    if ax2.is_first_row():
        ax2.set_xlabel( r'$z$', fontsize=24 )
        
save_file = 'r_vs_time_{}.pdf'.format( defaults['tag'] )
plotting.save_fig(
    out_dir = file_manager.get_project_figure_dir(),
    save_file = save_file,
    fig = fig,
)

fig


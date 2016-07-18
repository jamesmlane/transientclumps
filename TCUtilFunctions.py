import os
import sys
import numpy as np
import pdb
import subprocess
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from astropy.coordinates import Angle
from astropy import units

##### Plotting and Analysis Functions #####




### peakr_scatterplot ###
# Take the flux of an object and it's observered peak ratios to plot variation
# as a function of peak brightness. flux and peakr should contain all the
# measurements for all sources concatenated into one array each.
#
# Keywords:
# flux - Peak brightnesses for all sources (float array) -- required.
# peakr - calculated peak ratios for all sources (float array) -- required.
# y_hi - upper plotting range for y, should be calculated beforehand to prevent
# 		aliasing of bright or faint signals (float).
# y_lo - lower plotting range for y, same considerations as above (float).
# proto - protostellar classes for all objects (string array) -- required.
# solid - array denoting whether sources are solid (1, found in all maps) or
#		weak (0, not found in all maps) (int array) -- requred.
# orion - Are the sources in Orion A/B (boolean) -- required.

################################################################################
def peakr_scatterplot(flux, peakr, y_hi=1.25, y_lo=0.8, proto, solid, orion):

	#Set the x-range.
	x_hi = max(flux)*1.1
	x_lo = 0

	#Plot limits and labels.
	plt.xlabel('Peak Flux (Jy/Beam)')
	plt.xlim(x_lo,x_hi)
	plt.ylim(y_lo,y_hi)

	#Call the legend/scatterplot making function
	legend_scatter(flux, peakr, proto, solid, orion)

	#Plot reference lines.
	plt.plot((x_lo,x_hi),(1.0,1.0), '-k')
	plt.plot((x_lo,x_hi),(0.9,0.9), '--k')
	plt.plot((x_lo,x_hi),(1.1,1.1), '--k')
#def
################################################################################




### off_scatterplot ###
# Take the flux of an object and it's calculated peak offsets to plot
# as a function of peak brightness. flux and off should contain all the
# measurements for all sources concatenated into one array each.
#
# Keywords:
# flux - Peak brightnesses for all sources (float array) -- required.
# off - Calculated peak offsets for all sources (float array) -- required.
# y_hi - Upper plotting range for y (float).
# y_lo - Lower plotting range for y (float).
# proto - Protostellar classes for all objects (string array) -- required.
# solid - Array denoting whether sources are solid (1, found in all maps) or
#		weak (0, not found in all maps) (int array) -- requred.
# orion - Are the sources in Orion A/B (boolean) -- required.

################################################################################
def off_scatterplot(flux, off, y_hi=7, y_lo=-7, proto, solid, orion):

	#Set the x-range.
	x_hi = max(flux)*1.1
	x_lo = 0

	#Plot limits and labels.
	plt.xlabel('Peak Flux (Jy/Beam)')
	plt.ylim(y_lo,y_hi)
	plt.xlim(x_lo,x_hi)

	#Call the legend/scatterplot making function
	legend_scatter(flux, off, proto, solid, orion)

	#Plot a reference line.
	plt.plot((x_lo,x_hi),(0.0,0.0), '-k')
#def
################################################################################




### off_scatterplot ###
# Plot a time series of peak brightness for a single source.
#
# Keywords:
# time - Observation time in days (zeroed at 1st obs) (float array) -- required.
# peakr - Calculated peak brightness ratios (float array) -- required.
# flux - Source peak flux for annotation (float) -- required.
# proto - Protostellar class for annotation (string) -- required.
# num - Source number in base catalog (int) -- required.
# t_init - Date of first observation (time = 0) (string) -- required.
# posx - RA of source in degrees (float) -- required.
# posy- Declination of source in degrees (float) -- required.

################################################################################
def time_scatterplot(time, peakr, flux, proto, num, t_init, posx, posy):

	#Convert the decimal degrees to HMS strings.
	ra_deg = Angle(str(posx)+'d')
	dec_deg = Angle(str(posy)+'d')
	ra_hms = ra_deg.to_string(unit=units.hour, sep=':', precision=2)
	dec_dms = dec_deg.to_string(unit=units.degree, sep=':', precision=2)

	#Set the x-range.
	x_hi = max(time)*1.1
	x_lo = -max(time)*0.1
	plt.xlim(x_lo, x_hi)

	#Set the y-range.
	y_hi = 1.25
	y_lo = 0.75
	if np.amax(peakr) > 1.25:
		y_hi = np.amax(peakr)*1.5
	##fi
	if np.amax(peakr) < 0.75:
		y_lo = min(peakr)*0.7
	##fi
	plt.ylim(y_lo, y_hi)

	#Define the labels.
	plt.xlabel('Time (days + '+t_init+')')

	#Write important details about the source on the plot.
	plt.annotate('Source Number: '+str(num), xy = (0.05,0.2),
				 xycoords = 'axes fraction', size = 'small')
	plt.annotate('Peak Flux: '+'{:.4}'.format(flux)+' (Jy/Beam)',
				 xy = (0.05,0.15), xycoords = 'axes fraction', size = 'small')
	plt.annotate('Class: '+proto, xy = (0.05,0.1), xycoords = 'axes fraction',
				 size = 'small')
	plt.annotate('RA: '+ra_hms+',  Dec: '+dec_dms, xy = (0.05,0.05),
				 xycoords = 'axes fraction', size = 'small')

	#Plot the data.
	plt.scatter(time, peakr, alpha = 0.6, s = 60, marker = 'o', facecolor = 'r')

	#Plot reference lines.
	plt.plot((x_lo,x_hi),(1.0,1.0), '-k')
	plt.plot((x_lo,x_hi),(0.9,0.9), '--k')
	plt.plot((x_lo,x_hi),(1.1,1.1), '--k')
#def
################################################################################




### get_obs_date ###
# Returns the elapsed time in days for each observation in .sdf format.
#
# Keywords:
# img_name - Array of image names (in .sdf format) corresponding to images
# 			that will have the time calculated -- required.
#
# Returns:
# obs_time - Elapsed time between each observation and the first observation
# 			in days (float array).
# t_initial - The date of the first observation, corresponds with obs_time[0]
#			and is the zero point for obs_time[i] (string).
################################################################################
def get_obs_days(img_name):

	#If the user needs kappa called then do that.
	if call_kappa == True:
		subprocess.call('{ . $KAPPA_DIR/kappa.sh ;}', shell=True)
	##fi

	HST_end = np.empty(n_sdf, dtype='O')

	#Pipe the output of the fitsval call into an array.
	for i in range( len( img_name ) ):
		date_command = '$KAPPA_DIR/fitsval ndf='+img_name[i]+' keyword=HSTEND'
		cur_process = subprocess.Popen(date_command, shell=True,
										stdout=subprocess.PIPE)
		cur_stdout = cur_process.communicate()
		HST_end[i] = str(cur_stdout[0])[2:-3]

	obs_time = np.zeros(n_sdf)
	t_format = '%Y-%m-%dT%H:%M:%S'

	#Iterate over all observations, comparing their time to the first
	# observations to get an elapsed time in days.
	for i in range(n_sdf-1):
		delta_t = dt.strptime(HST_end[i+1], t_format) - dt.strptime(HST_end[0],
																	 t_format)
		#determine the total time in days
		days = delta_t.days
		seconds = delta_t.seconds/864000.0
		obs_time[i+1] = seconds+days
	###i

	#Format the starting date for further plotting.
	t_initial = dt.strptime(HST_end[0], t_format).strftime('%Y-%b-%d')

	return obs_time, t_initial

#def
################################################################################



### legend_scatter ###
# Plot sources and create a legend with colouring corresponding to the
# protostellar classification and shapes corresponding to the detection type.
#
# Keywords:
# flux - Peak flux of the sources being plotted (float array) -- required.
# variation - The variation being plotted (float array) -- required.
# proto - The protostellar classes of the sources (string array) -- required.
# solid - Array denoting whether sources are solid (1, found in all maps) or
#		weak (0, not found in all maps) (int array) -- requred.
# orion - Are the sources in Orion A/B (boolean) -- required.

################################################################################
def legend_scatter(flux, variation, proto, solid, orion):

	#If we are in Orion use Megeath nomenclature.
	if orion == 'TRUE':

		#Declare some tracking arrays.
		where_classP = np.where(proto == 'P')[0]
		where_classD = np.where(proto == 'D')[0]
		where_classFP = np.where(proto == 'FP')[0]
		where_classRP = np.where(proto == 'RP')[0]
		where_classNA = np.where(proto == 'NA')[0]
		where_solid = np.where(solid == 1)[0]
		where_weak = np.where(solid == 0)[0]
		legend_handle_check = np.zeros(7)

		#Check for protostars.
		if len(where_classP) > 0:
			#Figure out which class P sources are solid and weak.
			where_classP_solid = np.where(proto[where_solid] == 'P')[0]
			where_classP_weak = np.where(proto[where_weak] == 'P')[0]

			#Plot the solid and weak points.
			plt.scatter(flux[where_solid[where_classP_solid]],
						variation[where_solid[where_classP_solid]],
							 alpha=0.6, s=60, marker='o', facecolor='r')
			plt.scatter(flux[where_weak[where_classP_weak]],
						variation[where_weak[where_classP_weak]],
							 alpha=0.6, s=60, marker='^', facecolor='r')
			legend_handle_check[2] = 1

		##fi

		#Check for disks.
		if len(where_classD) > 0:
			#Figure out which class D sources are solid and weak.
			where_classD_solid = np.where(proto[where_solid] == 'D')[0]
			where_classD_weak = np.where(proto[where_weak] == 'D')[0]

			#Plot the solid and weak points.
			plt.scatter(flux[where_solid[where_classD_solid]],
					    variation[where_solid[where_classD_solid]],
							 alpha=0.6, s=60, marker='o', facecolor='b')
			plt.scatter(flux[where_weak[where_classD_weak]],
							 variation[where_weak[where_classD_weak]],
							 alpha=0.6, s=60, marker='^', facecolor='b')
			legend_handle_check[3] = 1
		##fi

		#Check for faint protostars.
		if len(where_classFP) > 0:
			#Figure out which class FP sources are solid and weak.
			where_classFP_solid = np.where(proto[where_solid] == 'FP')[0]
			where_classFP_weak = np.where(proto[where_weak] == 'F{P')[0]

			#Plot the solid and weak points.
			plt.scatter(flux[where_solid[where_classFP_solid]],
						variation[where_solid[where_classFP_solid]],
							 alpha=0.6, s=60, marker='o', facecolor='m')
			plt.scatter(flux[where_weak[where_classFP_weak]],
						variation[where_weak[where_classFP_weak]],
							 alpha=0.6, s=60, marker='^', facecolor='m')

			legend_handle_check[4] = 1
		##fi

		#Check for red protostars.
		if len(where_classRP) > 0:
			#Figure out which class RP sources are solid and weak.
			where_classRP_solid = np.where(proto[where_solid] == 'RP')[0]
			where_classRP_weak = np.where(proto[where_weak] == 'RP')[0]

			#Plot the solid and weak points.
			plt.scatter(flux[where_solid[where_classRP_solid]],
						variation[where_solid[where_classRP_solid]],
							 alpha=0.6, s=60, marker='o', facecolor='y')
			plt.scatter(flux[where_weak[where_classRP_weak]],
						variation[where_weak[where_classRP_weak]],
							 alpha=0.6, s=60, marker='^', facecolor='y')
			legend_handle_check[5] = 1
		##fi

		#Check for sources that aren't matched.
		if len(where_classNA) > 0:
			#Figure out which class NA sources are solid and weak
			where_classNA_solid = np.where(proto[where_solid] == 'NA')[0]
			where_classNA_weak = np.where(proto[where_weak] == 'NA')[0]

			#Plot the solid and weak points.
			plt.scatter(flux[where_solid[where_classNA_solid]],
						variation[where_solid[where_classNA_solid]],
							 alpha=0.6, s=60, marker='o', facecolor='none')
			plt.scatter(flux[where_weak[where_classNA_weak]],
						variation[where_weak[where_classNA_weak]],
							 alpha=0.6, s=60, marker='^', facecolor='none')
			legend_handle_check[6] = 1
		##fi

		#Check for solid sources (full identifications).
		if len(where_solid) > 0:
			#Scatter plot these sources
			plt.scatter(flux[where_solid], variation[where_solid], s=60,
						facecolor='none', edgecolor='k', marker='o')
			legend_handle_check[0] = 1
		##fi

		#Check for weak sources (partial identifications).
		if len(where_weak) > 0:
			plt.scatter(flux[where_weak], variation[where_weak], s=60,
						facecolor='none', edgecolor='k', marker='^')
			legend_handle_check[1] = 1
		##fi

		#Make all of the legend handles.
		black_circle = mlines.Line2D((0,0),(0,0), marker = 'o',
									 mec = 'k', mfc = 'none', mew = 1.2,
									 ms = 10, linestyle = '', label='Solid')
		black_triangle = mlines.Line2D((0,1),(0,0),  marker = '^',
									   mec = 'k', mfc = 'none', mew = 1.2,
									   ms = 10, linestyle = '', label='Partial')
		red_patch = mpatches.Patch(fc= 'r', ec = 'k', label='P')
		blue_patch = mpatches.Patch(fc = 'b', ec = 'k', label='D')
		magenta_patch = mpatches.Patch(fc = 'm', ec = 'k', label='FP')
		yellow_patch = mpatches.Patch(fc = 'y', ec = 'k', label='RP')
		white_patch = mpatches.Patch(fc = 'none', ec = 'k', label='N/A')

		#Make an array of objects to hold the legend handles.
		legend_handle_temp = np.array([black_circle, black_triangle, red_patch,
										blue_patch, magenta_patch, yellow_patch,
										white_patch],dtype='O')

		#Choose the legend handles that correspond to the identified classes.
		legend_handle_actual = legend_handle_temp[np.where(legend_handle_check
														   == 1)[0]]

		#Return the current axis object to change the bounds.
		ax = plt.gca()
		cur_bounds = ax.get_position()
		ax.set_position([cur_bounds.x0, cur_bounds.y0, cur_bounds.width * 0.8,
						 cur_bounds.height])

		#Place the legend off to one side.
		plt.legend(handles=list(legend_handle_actual), bbox_to_anchor=(1, 0.5),
				   loc='center left')
	##fi

	#If we aren't in Orion then use Dunham nomenclature.
	if orion == 'FALSE':

		#Declare some tracking arrays.
		where_class01 = np.where(proto == '0+1')[0]
		where_classF = np.where(proto == 'F')[0]
		where_class2 = np.where(proto == '2')[0]
		where_class3 = np.where(proto == '3')[0]
		where_classNA = np.where(proto == 'NA')[0]
		where_solid = np.where(solid == 1)[0]
		where_weak = np.where(solid == 0)[0]
		legend_handle_check = np.zeros(7)

		#Check for class 0+1 sources.
		if len(where_class01) > 0:
			#Figure out which class 01 sources are solid and weak.
			where_class01_solid = np.where(proto[where_solid] == '0+1')[0]
			where_class01_weak = np.where(proto[where_weak] == '0+1')[0]

			#Plot the solid and weak points.
			plt.scatter(flux[where_solid[where_class01_solid]],
						variation[where_solid[where_class01_solid]],
							 alpha=0.6, s=60, marker='o', facecolor='r')
			plt.scatter(flux[where_weak[where_class01_weak]],
						variation[where_weak[where_class01_weak]],
							 alpha=0.6, s=60, marker='^', facecolor='r')
			legend_handle_check[2] = 1
		##fi

		#Check for flat spectrum sources.
		if len(where_classF) > 0:
			#Figure out which class F sources are solid and weak.
			where_classF_solid = np.where(proto[where_solid] == 'F')[0]
			where_classF_weak = np.where(proto[where_weak] == 'F')[0]

			#Plot the solid and weak points.
			plt.scatter(flux[where_solid[where_classF_solid]],
					    variation[where_solid[where_classF_solid]],
							 alpha=0.6, s=60, marker='o', facecolor='b')
			plt.scatter(flux[where_weak[where_classF_weak]],
							 variation[where_weak[where_classF_weak]],
							 alpha=0.6, s=60, marker='^', facecolor='b')
			legend_handle_check[3] = 1
		##fi

		#Check for class 2 sources
		if len(where_class2) > 0:
			#Figure out which class 2 sources are solid and weak.
			where_class2_solid = np.where(proto[where_solid] == '2')[0]
			where_class2_weak = np.where(proto[where_weak] == '2')[0]

			#Plot the solid and weak points.
			plt.scatter(flux[where_solid[where_class2_solid]],
						variation[where_solid[where_class2_solid]],
							 alpha=0.6, s=60, marker='o', facecolor='m')
			plt.scatter(flux[where_weak[where_class2_weak]],
						variation[where_weak[where_class2_weak]],
							 alpha=0.6, s=60, marker='^', facecolor='m')
			legend_handle_check[4] = 1
		##fi

		#Check for class 3 sources.
		if len(where_class3) > 0:
			#Figure out which class 3 sources are solid and weak.
			where_class3_solid = np.where(proto[where_solid] == '3')[0]
			where_class3_weak = np.where(proto[where_weak] == '3')[0]

			#Plot the solid and weak points.
			plt.scatter(flux[where_solid[where_class3_solid]],
						variation[where_solid[where_class3_solid]],
							 alpha=0.6, s=60, marker='o', facecolor='y')
			plt.scatter(flux[where_weak[where_class3_weak]],
						variation[where_weak[where_class3_weak]],
							 alpha=0.6, s=60, marker='^', facecolor='y')
			legend_handle_check[5] = 1
		##fi

		#Check for class NA sources.
		if len(where_classNA) > 0:
			#Figure out which unmatched sources are solid and weak.
			where_classNA_solid = np.where(proto[where_solid] == 'NA')[0]
			where_classNA_weak = np.where(proto[where_weak] == 'NA')[0]

			#Plot the solid and weak points.
			plt.scatter(flux[where_solid[where_classNA_solid]],
						variation[where_solid[where_classNA_solid]],
							 alpha=0.6, s=60, marker='o', facecolor='none')
			plt.scatter(flux[where_weak[where_classNA_weak]],
						variation[where_weak[where_classNA_weak]],
							 alpha=0.6, s=60, marker='^', facecolor='none')
			legend_handle_check[6] = 1
		##fi

		#Check for solid sources (full identifications).
		if len(where_solid) > 0:
			#Scatter plot these sources.
			plt.scatter(flux[where_solid], variation[where_solid], s=60,
						color='none', edgecolor='k', marker='o')
			legend_handle_check[0] = 1
		##fi

		#Check for weak sources (partial identifications).
		if len(where_weak) > 0:
			#Scatter plot these sources.
			plt.scatter(flux[where_weak], variation[where_weak], s=60,
						color='none', edgecolor='k', marker='^')
			legend_handle_check[1] = 1
		##fi

		#Make all of the legend handles.
		black_circle = mlines.Line2D((0,0),(0,0), marker = 'o',
									 mec = 'k', mfc = 'none', mew = 1.2,
									 ms = 10, linestyle = '', label='Solid')
		black_triangle = mlines.Line2D((0,1),(0,0),  marker = '^',
									   mec = 'k', mfc = 'none', mew = 1.2,
									   ms = 10, linestyle = '', label='Partial')
		red_patch = mpatches.Patch(fc = 'r', ec = 'k', label='0+1')
		blue_patch = mpatches.Patch(fc = 'b', ec = 'k', label='F')
		magenta_patch = mpatches.Patch(fc = 'm', ec = 'k', label='2')
		yellow_patch = mpatches.Patch(fc = 'y', ec = 'k', label='3')
		white_patch = mpatches.Patch(fc = 'none', ec = 'k', label='N/A')

		#Make an array of objects to hold the legend handles.
		legend_handle_temp = np.array([black_circle, black_triangle, red_patch,
										blue_patch, magenta_patch, yellow_patch,
										white_patch],dtype='O')

		#Choose the legend handles that correspond to the identified classes.
		legend_handle_actual = legend_handle_temp[np.where(legend_handle_check
														   == 1)[0]]

		#Return the current axis object to change the bounds.
		ax = plt.gca()
		cur_bounds = ax.get_position()
		ax.set_position([cur_bounds.x0, cur_bounds.y0, cur_bounds.width * 0.8,
						 cur_bounds.height])

		#Place the legend off to one side.
		plt.legend(handles=list(legend_handle_actual), bbox_to_anchor=(1, 0.5),
				   loc='center left')
	##fi
#def
################################################################################

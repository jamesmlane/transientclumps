import numpy as np
import scipy.spatial.distance as ssd
import astropy.io.fits as apfits

### source_match ###
# Takes a gaussclumps catalog for a target image and compares the catalog
# to a reference catalog to find matches and report offsets in position and
# flux.
#
# Keywords:
# cat_name - name of the target gaussclumps catalog in .fits format (string).
#          -- required
# ref_name - name of the reference gaussclumps catalog in .fits format (string).
#          -- required
# minpeak - minimum brightness of a source to consider matching it (Jy/Beam).
# maxrad - maximum effective radius of a source to consider matching it (").
# maxsep - maximum allowed separation between two sources to match them (").
#
# Output:
# off - a 4-element array that gives the x and y offset of the target image
#           with respect to the reference image (ie: target - reference) as well
#           as peak and radius ratios (ie: target / reference). All values are
#           averages.
# err - a 4-element array that gives the standard deviation in each of the above
#           measurements, respectively.

def source_match(cat_name, ref_name, minpeak=0.2, maxrad=30, maxsep=10):

    #Read in the reference catalog and prepare the data.
    ref_data = apfits.getdata(ref_name, 0)
    ref_nsources = len(ref_data['Cen1'])
    ref_posxvals = ref_data['Cen1']
    ref_posyvals = ref_data['Cen2']
    ref_peakvals = ref_data['Peak']
    ref_fwhm1vals = ref_data['GCFWHM1']
    ref_fwhm2vals = ref_data['GCFWHM2']

    #Calculate the effective radius, the square root of the two FWHM multiplied
    # together. Multiply by 1.5 to account for FWHM is diameter but want radius,
    # and units are pixels but want arcseconds.
    ref_rvals = np.sqrt( np.multiply( ref_fwhm1vals, ref_fwhm2vals ) ) * 1.5

    #Keep an array to track index of successful matches in the target catalog.
    matched_sources = np.zeros(ref_nsources)
    matched_sources[:] = -1

    #Read in the target catalog and prepare the data.
    targ_data = apfits.getdata(cat_name, 0)
    targ_nsources = len(targ_data['Cen1'])
    targ_posxvals = targ_data['Cen1']
    targ_posyvals = targ_data['Cen2']
    targ_peakvals = targ_data['Peak']
    targ_fwhm1vals = targ_data['GCFWHM1']
    targ_fwhm2vals = targ_data['GCFWHM2']

    #Calculate effective radius, same as above.
    targ_rvals = np.sqrt( np.multiply( targ_fwhm1vals, targ_fwhm2vals ) ) * 1.5

    #Keep an array to track whether or not target catalog sources have found
    # a match, because they can't match more than once. Brightest sources
    # get priority.
    targ_matchcheck = np.zeros(targ_nsources)

    #Loop over all sources in the reference catalog.
    for i in range(ref_nsources):

        #Test a radius and peak condition here:
        if ref_rvals[i] > maxrad: continue
        if ref_peakvals[i] < minpeak: continue

        #Declare the reference coordinates.
        ref_coords = np.array([[ref_posxvals[i], ref_posyvals[i]],])

        #Set a minimum distance tracking variable and the matched source
        # indexing variable.
        mindist = -1
        match_ind = -1

        #Loop over the target catalog sources.
        for j in range(targ_nsources):

            #Don't match a source that already has a match.
            if targ_matchcheck[j] == 1: continue

            #Declare the target coordinates.
            targ_coords = np.array([[targ_posxvals[j], targ_posyvals[j]],])

            #Calculate distance.
            dist = ssd.cdist(ref_coords, targ_coords, metric='euclidean')

            #Check if this source meets the criterion for matching distance.
            if dist < maxsep:

                #If this is the first iteration then mindist will be -1
                if mindist == -1:
                    #Auto assign mindist and match_ind
                    match_ind = j
                    mindist = dist
                ###fi

                #Now check if a new closest source.
                if dist < mindist:
                    #Remember the new source number and assign dist to mindist.
                    match_ind = j
                    mindist = dist
                    ###fi
            ###fi
        ###j

        #If we are at the end of the reference catalog then assign
        # this source if it found a match.
        if match_ind >= 0:
            targ_matchcheck[match_ind] = 1
            matched_sources[i] = match_ind
        ##fi
    ###i

    #Find out the number of matched sources.
    n_matched = len(np.where(matched_sources != -1)[0])

    #Calculate averages, ratios, and offsets. Only if sources matched.
    if n_matched > 0:
        matched_peak_ref = np.zeros(n_matched)
        matched_posx_ref = np.zeros(n_matched)
        matched_posy_ref = np.zeros(n_matched)
        matched_r_ref = np.zeros(n_matched)

        #Make ratio and offset catalogs.
        matched_peakr = np.zeros(n_matched)
        matched_xoff = np.zeros(n_matched)
        matched_yoff = np.zeros(n_matched)
        matched_rr = np.zeros(n_matched)
        match_count = 0
        matched_index = np.zeros((n_matched,2), dtype='int')

        #Loop over the number of sources in the reference catalog.
        for i in range(ref_nsources):
            #Decide whether or not this source was matched or not.
            if matched_sources[i] == -1: continue
            ##fi
            #Then it's a good source.
            targ_match = int(matched_sources[i])
            matched_peak_ref[match_count] = ref_peakvals[i]
            matched_posx_ref[match_count] = ref_posxvals[i]
            matched_posy_ref[match_count] = ref_posyvals[i]
            matched_r_ref[match_count] = ref_rvals[i]

            #Create offsets and ratios (multiple by 3600 for ")
            matched_xoff[match_count] = (targ_posxvals[targ_match]-
            ref_posxvals[i])*3600.0
            matched_yoff[match_count] = (targ_posyvals[targ_match]-
            ref_posyvals[i])*3600.0
            matched_peakr[match_count] = (targ_peakvals[targ_match]/
            ref_peakvals[i])
            matched_rr[match_count] = (targ_rvals[targ_match]/
            ref_rvals[i])

            match_count+=1
        ###i

        #Calculate averages.
        avg_xoff = np.average(matched_xoff)
        avg_yoff = np.average(matched_yoff)
        avg_peakr = np.average(matched_peakr)
        avg_rr = np.average(matched_rr)

        #Calculate error.
        std_xoff = np.std(matched_xoff)
        std_yoff = np.std(matched_yoff)
        std_peakr = np.std(matched_peakr)
        std_rr = np.std(matched_rr)

        print('\n'+str(match_count)+' sources matched.\n')
        return [avg_xoff,avg_yoff,avg_peakr,avg_rr], [std_xoff,std_yoff,std_peakr,std_rr]
    ##fi

    #If no sources found return None.
    if n_matched == 0:
        print('\nNo sources found.\n')
        return [None,None,None,None], [None,None,None,None]
#def

"""
 ------------------
 READ_SPICAV_IR.PY

 This file contains the tools necessary to read and store the SPICAV IR data.

 Based on Anna Fedorova's Matlab version. Translated into Python,
 corrected and completed by Loic Rossi. 2013, 2014, 2015
 Authors : Loic Rossi, Gourav Mahapatra
 Licence : GNU/GPL and CeCILL
 2012 - 2015

 CONTENT
   Definitions of classes:
      Geo
      Data
   Functions:
      read_radiance(filename)
      read_geometry(filename_geo)
      adjust_radiance0(w0,w1,r0,r1)
      get_bande
      mean_data_phase
      clean_nan(array_in)
      processing_orbits(orbit_n,orbit_a,...)
"""
# --- Import modules ---

#import struct as st
import numpy as np
import matplotlib.pyplot as mpl
#from fovint import fovint
#import array
import os

#from decimal import Decimal

# ------ CLASS DEFINITIONS --------

class Geo:
    """ This class defines an object containing the geometric data of a
    SPICAV-IR observation

    Attributes:
            date: date of observation
            time: timestamp of each point
            time_str: timestamp string of each point
            L_sun: solar longitude
            Lat_sc: spacecraft latitude
            Lon_sc: spacecraft longitude
            Lat: Latitude of obs point on planet
            Lon: longitude
            Alt_sc: altitude of spacecraft
            H: altitude of the closest point of the LOS to the surface. Used
                for stellar occultations. Should be negative in nadir (LOS goes
                in the planet).
            phase: phase angle
            sza: solar zenith angle
            emission: emission angle
            local_time: local time of obs pt
            S_Mdist: distance btw Sun and Venus
            target_dist: distance to the targeted point on the planet
            alt_topo: topographic altitude

    Methods:
        __init__
        __repr__
        __add__
        calc_azimuth

    """

    def __init__(self, date = [], time = [], time_str = [], L_sun = [],
                 Lat_sc = [], Lon_sc = [], Lat = [], Lon = [] ,
                 Alt_sc = [],H = [], phase = [], sza = [], emission = [],
                 local_time = [], S_Mdist = [],
                 target_dist = [], alt_topo = []):

        """ Initialisation of the geometric data from SPICAV-IR observation."""

        self.date = date
        self.time = time
        self.time_str = time_str
        self.L_sun = L_sun
        self.Lat_sc = Lat_sc
        self.Lon_sc = Lon_sc
        self.Lat = Lat
        self.Lon = Lon
        self.Alt_sc = Alt_sc
        self.H = H
        self.phase = phase
        self.sza = sza
        self.emission = emission
        self.local_time = local_time
        self.S_Mdist = S_Mdist
        self.target_dist = target_dist
        self.alt_topo = alt_topo

    def __repr__(self):
        """Displays some basic geometric information about the observation."""

        str = " Phase coverage: {} to {}".format(min(self.phase),
                                                 max(self.phase))
        str = str + "\n Latitude: {} to {}".format(min(self.Lat),
                                                   max(self.Lat))
        str = str + "\n local time: {} to {}".format(min(self.local_time),
                                                     max(self.local_time))

        return str

    def __add__(self,geom):
        """
        This method defines the addition of two geometrical objects. The
        resulting object will contain informations from both initial objects.

        Example:
            new_geo = geo1 + geo2
        """
        new_geo = Geo()

        if self.phase==[]:  # check if geometric object is blank
            new_geo.date = geom.date
            new_geo.time = geom.time
            new_geo.time_str = geom.time_str
            new_geo.L_sun = geom.L_sun
            new_geo.Lat_sc = geom.Lat_sc
            new_geo.Lon_sc = geom.Lon_sc
            new_geo.Lat = geom.Lat
            new_geo.Lon = geom.Lon
            new_geo.Alt_sc = geom.Alt_sc
            new_geo.H = geom.H
            new_geo.phase = geom.phase
            new_geo.sza = geom.sza
            new_geo.emission = geom.emission
            new_geo.local_time = geom.local_time
            new_geo.S_Mdist = geom.S_Mdist
            new_geo.target_dist = geom.target_dist
            new_geo.alt_topo = geom.alt_topo
        else:
            new_geo.date = np.append(self.date,geom.date,axis=1)
            new_geo.time = np.append(self.time,geom.time,axis=1)
            new_geo.time_str = np.append(self.time_str,geom.time_str,axis=0)
            new_geo.L_sun = np.append(self.L_sun,geom.L_sun)
            new_geo.Lat_sc = np.append(self.Lat_sc,geom.Lat_sc)
            new_geo.Lon_sc = np.append(self.Lon_sc,geom.Lon_sc)
            new_geo.Lat = np.append(self.Lat,geom.Lat)
            new_geo.Lon = np.append(self.Lon,geom.Lon)
            new_geo.Alt_sc = np.append(self.Alt_sc,geom.Alt_sc)
            new_geo.H = np.append(self.H,geom.H)
            new_geo.phase = np.append(self.phase,geom.phase)
            new_geo.sza = np.append(self.sza,geom.sza)
            new_geo.emission = np.append(self.emission,geom.emission)
            new_geo.local_time = np.append(self.local_time,geom.local_time)
            new_geo.S_Mdist = np.append(self.S_Mdist,geom.S_Mdist)
            new_geo.target_dist = np.append(self.target_dist,geom.target_dist)
            new_geo.alt_topo = np.append(self.alt_topo,geom.alt_topo)

        return new_geo


    def calc_azimuth(self):
        """
        This method computes the azimuth angle from the geometric
        data. To be used once all geo data have been read and treated.

        Usage:
            data.geo.calc_azimuth()

        """
        theta = np.radians(self.sza)
        thetap = np.radians(self.emission)
        alpha = np.radians(self.phase)
        c_delta_phi = ( (np.cos(alpha) - (np.cos(theta)*np.cos(thetap)) )
                       /(np.sin(theta)*np.sin(thetap)) )

        delta_phi = np.degrees(np.arccos(c_delta_phi))

        self.azimuth = delta_phi


class Data:
    """ This class contains the main data from an observation: radiances,
    wavelengths, polarization...
    includes the geometric informations, as a sub class

    Attributes:
        w0: wavelengths for detector 0
        w1: wavelengths for detector 1
        r0: radiances for detector 0
        r1: radiances for detector 1
        wd0: averaged DOTS wavelengths for detector 0
        wd1: averaged DOTS wavelengths for detector 1
        rd0: averaged DOTS radiances for detector 0
        rd1: averaged DOTS radiances for detector 1
        time: instant of the beginning of observation in format
            [YYYY, MM, DD, hh, mm, ss, .s]
        pol: array containing the polarisation of the DOTS.
            Shape (number of DOTS, number of measures)
        dpol: uncertainty of measure on pol. same shape
        intensity: radiance (W/m2/sr/um) measured for DOTS. shaped as pol
        dI: uncertainty on the intensity
        geo: a Geo object
        phase: array with phase angle. Identical to geo.phase. Here to have an
            easier access to the phase angle.
        num_orbit: array with orbit number and observation sequence.
            format [[ORB, SEQ], [ORB,SEQ],...] shape (nb of obs, 2)
    Methods:
        __init__
        adjust_radiance
        filter
        split
        __iadd__
        __add__

    """

    def __init__(self,w0 = [],w1 = [],r0 = [],r1 = [],wd0 = [],wd1 = []
                 ,rd0 = [],rd1 = [],time = [],pol = [],dpol = [],
                 intensity = [], dI = [], geo = [],phase = [],num_orbit = []):
        """Initialisation of the data object."""
        self.w0 = w0
        self.w1 = w1
        self.r0 = r0
        self.r1 = r1
        self.wd0 = wd0
        self.wd1 = wd1
        self.rd0 = rd0
        self.rd1 = rd1
        self.time = time
        self.pol = pol
        self.dpol = dpol
        self.intensity = intensity
        self.dI = dI
        self.geo = geo
        self.phase = phase
        self.num_orbit= num_orbit

    def __repr__(self):
        """Displays some basic information about the data object"""

        if len(np.shape(self.time))>1:
            num_orbit_i = self.num_orbit[0].astype(int)
            time_i = self.time[0].astype(int)

            num_orbit_e = self.num_orbit[-1].astype(int)
            time_e = self.time[-1].astype(int)

            str = ('This data set is from \n'
                   'Orbit '
                   '{num_orbit[0]}N'
                   '{num_orbit[1]}\n').format(num_orbit=num_orbit_i)
            str = str + ('Date:'
                         '{time[0]}-{time[1]}-{time[2]}'
                         'T{time[3]}:{time[4]}:{time[5]}.{time[6]}'
                         '\n').format(time=time_i)

            str = str + ('to Orbit '
                         '{num_orbit[0]}N{num_orbit[1]}'
                         '\n').format(num_orbit=num_orbit_e)
            str = str + ('Date: {time[0]}-{time[1]}-{time[2]}'
                         'T{time[3]}:{time[4]}:{time[5]}.{time[6]}'
                         '\n').format(time=time_e)
        else:
            num_orbit_i = self.num_orbit[0].astype(int)
            time_i = self.time.astype(int)

            str = ('This data set is from\n'
                   'Orbit {num_orbit[0]}N{num_orbit[1]}'
                   '\n').format(num_orbit=num_orbit_i)
            str = str + ('Date: {time[0]}-{time[1]}-{time[2]}'
                         'T{time[3]}:{time[4]}:{time[5]}.{time[6]}'
                         ' ').format(time=time_i)


        return str

    def adjust_radiance(self,lw_noise='lw_NEB_g2_t28.txt',
                        sw_noise='sw_NEB_g2_t28.txt'):

        """ ADJUST_RADIANCE

        This method adjusts the data arrays. It reduces the number of points in
        SPICAV IR observations to a given set of reference wavelengths.  It
        produces an averaged array of wavelenghts of shape (nb of ref. wvls, nb
        of obs) and produces the associated radiance array, and the
        polarization array obtained from the two detectors.

        INPUTS :
        --- the data object
        --- lw_noise: path to the file describing the noise equivalent
            brightness of the detectors in the LW band
        --- sw_noise: path to the file describing the noise equivalent
            brightness of the detectors in the SW band

        OUTPUTS : The method will generate new attributes :
        --- wd0 : array containing the wavelengths averaged by groups
        --- wd1 : array containing the wavelengths averaged by groups
        --- rd0 : array containing the radiance averaged by groups
        --- rd1 : array containing the radiance averaged by groups
        --- pol : polarization degree for each reference wavelength and for
            each observation
        --- dpol : the error on the polarization

        This function uses the reference frequencies of the dots mode of SPICAV
        to calculate the wavelengths of the dots.  The frequency step between
        two points in one dot is 32kHz.  Averages the data within reference
        dots (over 3 or 5 or 10 points / wvl) for all observation points.
        """

        #Initialisation
        wl_d0 = self.w0
        wl_d1 = self.w1
        radf_d0 = self.r0/2.
        radf_d1 = self.r1/2.
        # necessary to divide by two because of calibration used by Anna

        nobs = len(wl_d0[0,:])

        r0_invalid = (radf_d0<0)  # identification of negative radiance values
        r1_invalid = (radf_d1<0)

        radf_d0[r0_invalid] = np.nan  # negative values are set to NaN
        radf_d1[r1_invalid] = np.nan

        ref_freq_sw = np.array([241267,201654,
                                176169,150643])  # ref freq in kHz of the aotf SW

        ref_freq_lw = np.array([133414,126219,
                                122008,117377,
                                114260,111382,
                                109753,95137,92663,88359])

        #ref freq LW
        freq_step = 32. #step in kHz used for dots

        ref_freq_sw_p = ref_freq_sw + 5*freq_step
        # ref freq plus 5 freq increments
        ref_freq_lw_p = ref_freq_lw + 5*freq_step

        # Coefficients for SW and LW for the calibration polynomial. Line 0 is
        # det0, line 1 in det1.
        coef_sw = np.array([[-4.9405101e-8, 7.6969006e-2, -2.9822051e2],
                            [-5.0454785e-8, 7.7358519e-2, -3.3244465e2]])
        coef_lw = np.array([[-3.3865473e-8, 7.2595705e-2, -2.0449838],
                            [-3.5371703e-8, 7.2919764e-2, -1.9140569e1]])

        # Construction of the polynomial v = af**2 + bf + c
        nu_sw_0 = np.poly1d(coef_sw[0,:])  # function to get the wavenb for det0
        nu_sw_1 = np.poly1d(coef_sw[1,:])  # function to get the wavenb for det1
        nu_lw_0 = np.poly1d(coef_lw[0,:])
        nu_lw_1 = np.poly1d(coef_lw[1,:])

        # Application of the polynomial to the reference freq for both channels
        # and both detectors
        ref_wvl_sw_0 = 1./nu_sw_0(ref_freq_sw)  # calc to get the wvl for det0
        ref_wvl_sw_1 = 1./nu_sw_1(ref_freq_sw)  # calc to get the wvl for det1
        ref_wvl_lw_0 = 1./nu_lw_0(ref_freq_lw)
        ref_wvl_lw_1 = 1./nu_lw_1(ref_freq_lw)

        # Application of the polynomial to the reference freq. shifted by 5
        # increments
        ref_wvl_sw_p0 = 1./nu_sw_0(ref_freq_sw_p)
            #calc to get wvl plus 5 increments in freq det0
        ref_wvl_sw_p1 = 1./nu_sw_1(ref_freq_sw_p)
            #calc to get wvl plus 5 increments in freq det1
        ref_wvl_lw_p0 = 1./nu_lw_0(ref_freq_lw_p)
        ref_wvl_lw_p1 = 1./nu_lw_1(ref_freq_lw_p)

        # conversion of wavenumber into wavelength
        ref_wvl_sw = 1e7 * np.vstack((ref_wvl_sw_0,ref_wvl_sw_1))
            #wavenumber in cm-1 so we need to convert to nanometers
        ref_wvl_lw = 1e7 * np.vstack((ref_wvl_lw_0,ref_wvl_lw_1))
        ref_wvl_sw_p = 1e7 * np.vstack((ref_wvl_sw_p0,ref_wvl_sw_p1))
        ref_wvl_lw_p = 1e7 * np.vstack((ref_wvl_lw_p0,ref_wvl_lw_p1))

        ref_wvl = np.hstack((ref_wvl_sw,ref_wvl_lw)) #all ref wvl together
        ref_wvl_p = np.hstack((ref_wvl_sw_p,ref_wvl_lw_p)) #same plus increments

        #Noise equivalent brightness
        neb_sw = np.genfromtxt(sw_noise,skip_header=1)
        neb_lw = np.genfromtxt(lw_noise,skip_header=1)

        neb = np.vstack((neb_lw,neb_sw))

        # creation of final arrays
        avg_wvl0 = np.zeros([14,nobs]) # shape(nb of wvl, number of measures)
        avg_wvl1 = np.zeros([14,nobs])

        avg_rad0 = np.zeros([14,nobs])
        avg_rad1 = np.zeros([14,nobs])

        neb_det0 = np.zeros([14])
        neb_det1 = np.zeros([14])
        dpol = np.zeros([14,nobs])
        dI = np.zeros([14,nobs])

        # Averaging the data by wavelength dots for all observation points
        for i,wvl in enumerate(ref_wvl[0,:]):  #loop on wavelenghts

            # selection of wavelengths in the NEB data
            noise_idx = (neb[:,0]<ref_wvl[0,i]) * (neb[:,0]>ref_wvl_p[0,i])
            neb_det0[i] = np.mean(neb[noise_idx,1])
                #mask on the NEB data so we get the noise on detectors for
            neb_det1[i] = np.mean(neb[noise_idx,2])  # all reference wavelenghts

            for j,z in enumerate(wl_d0[0,:]):  #loop on observation times

                # identification of wavelengths by using a mask (excluding nans
                # on radiance))
                indices0 = np.where( (wl_d0[:,j]>ref_wvl_p[0,i])
                                    & (wl_d0[:,j]<ref_wvl[0,i])
                                    & (np.isfinite(radf_d0[:,j])) )
                indices1 = np.where( (wl_d1[:,j]>ref_wvl_p[1,i])
                                    & (wl_d1[:,j]<ref_wvl[1,i])
                                    & (np.isfinite(radf_d1[:,j])) )

                avg_wvl0[i,j] = np.median( wl_d0[indices0,j] )
                avg_wvl1[i,j] = np.median( wl_d1[indices1,j] )

                avg_rad0[i,j] = np.median( radf_d0[indices0,j] )
                avg_rad1[i,j] = np.median( radf_d1[indices1,j] )

        # polarization is calculated as (det1 - det0)/(det0 + det1)
        pol = (avg_rad1 - avg_rad0) / (avg_rad0 + avg_rad1)

        #calculus of the error on polarization for all observation points
        for i,d in enumerate(avg_rad0[0,:]):
            dpol[:,i] = np.sqrt( ((2*avg_rad0[:,i]*neb_det1)
                                  /(avg_rad0[:,i]+avg_rad1[:,i])**2)**2
                                +( (2*avg_rad1[:,i]*neb_det0)
                                  /(avg_rad0[:,i] + avg_rad1[:,i])**2)**2 )

            dI[:,i] = np.sqrt(neb_det1**2 + neb_det0**2)

        intensity = (avg_rad0 + avg_rad1)

        self.wd0 = avg_wvl0
        self.wd1 = avg_wvl1
        self.rd0 = avg_rad0
        self.rd1 = avg_rad1
        self.pol = pol
        self.dpol = dpol
        self.intensity = intensity
        self.dI = dI

        return self

    def __iadd__(self,new_data):
        """Defines the incrementation of a data objects.  The object
        contains all the initial informations plus the content of the
        new data object.  Warning: the append mode is not completely
        fonctionnal. Some bug may occur.

        Usage:
            data += more_data
        """

        if self.pol==[]:
            self.pol = new_data.pol
            self.dpol = new_data.dpol
            self.intensity = new_data.intensity
            self.dI = new_data.dI
            self.phase = new_data.geo.phase
            self.geo = new_data.geo #also adds the geometric info
            self.time = new_data.time
            self.num_orbit = new_data.num_orbit
            self.w0 = new_data.w0
            self.w1 = new_data.w1
            self.wd0 = new_data.wd0
            self.wd1 = new_data.wd1
            self.r0 = new_data.r0
            self.r1 = new_data.r1
            self.rd0 = new_data.rd0
            self.rd1 = new_data.rd1
        else:
            self.pol = np.append(self.pol,new_data.pol,axis=1)
            self.dpol = np.append(self.dpol,new_data.dpol,axis=1)
            self.intensity = np.append(self.intensity,new_data.intensity,
                                       axis=1)
            self.dI = np.append(self.dI,new_data.dI, axis=1)
            self.phase = np.append(self.phase,new_data.geo.phase)
            self.geo = self.geo.__add__(new_data.geo) #also adds the geometric info
            self.time = np.vstack((self.time,new_data.time))
            self.num_orbit = np.vstack((self.num_orbit,new_data.num_orbit))
#           if np.shape(self.w0)[0] == np.shape(new_data.w0)[0]: # BUG HERE
#               # pb with different array shapes
#               self.w0 = np.hstack((self.w0,new_data.w0))
#               self.w1 = np.hstack((self.w1,new_data.w1))
#               self.wd0 = np.hstack((self.wd0,new_data.wd0))
#               self.wd1 = np.hstack((self.wd1,new_data.wd1))
#               self.r0 = np.hstack((self.r0,new_data.r0))
#               self.r1 = np.hstack((self.r1,new_data.r1))
#               self.rd0 = np.hstack((self.rd0,new_data.rd0))
#               self.rd1 = np.hstack((self.rd1,new_data.rd1))

        return self


    def __add__(self,new_data):
        """Defines the addition of two data objects.
        The new object contains all the informations contained in the initial
        data objects.
        Warning: the append mode is not completely fonctionnal. Some bug may
        occur.

        Usage:
            new_data = data1 + data2
        """
        out_data = Data()

        if self.pol==[]:
            out_data.pol = new_data.pol
            out_data.dpol = new_data.dpol
            out_data.intensity = new_data.intensity
            out_data.dI = new_data.dI
            out_data.phase = new_data.geo.phase
            out_data.geo = new_data.geo #also adds the geometric info
            out_data.time = new_data.time
            out_data.num_orbit = new_data.num_orbit
            out_data.w0 = new_data.w0
            out_data.w1 = new_data.w1
            out_data.wd0 = new_data.wd0
            out_data.wd1 = new_data.wd1
            out_data.r0 = new_data.r0
            out_data.r1 = new_data.r1
            out_data.rd0 = new_data.rd0
            out_data.rd1 = new_data.rd1
        else:
            out_data.pol = np.append(self.pol,new_data.pol,axis=1)
            out_data.dpol = np.append(self.dpol,new_data.dpol,axis=1)
            out_data.intensity = np.append(self.intensity,new_data.intensity,
                                           axis=1)
            out_data.dI = np.append(self.dI,new_data.dI,
                                           axis=1)
            out_data.phase = np.append(self.phase,new_data.geo.phase)
            out_data.geo = self.geo.__add__(new_data.geo) #also adds the geometric info
            out_data.time = np.vstack((self.time,new_data.time))
            out_data.num_orbit = np.vstack((self.num_orbit,new_data.num_orbit))
            if np.shape(self.w0)[0] == np.shape(new_data.w0)[0]: # BUG HERE
                # pb with different array shapes
                out_data.w0 = np.hstack((self.w0,new_data.w0))
                out_data.w1 = np.hstack((self.w1,new_data.w1))
                out_data.wd0 = np.hstack((self.wd0,new_data.wd0))
                out_data.wd1 = np.hstack((self.wd1,new_data.wd1))
                out_data.r0 = np.hstack((self.r0,new_data.r0))
                out_data.r1 = np.hstack((self.r1,new_data.r1))
                out_data.rd0 = np.hstack((self.rd0,new_data.rd0))
                out_data.rd1 = np.hstack((self.rd1,new_data.rd1))

        return out_data

    def add_geometry(self,filename_geo):
        """ Creates a Geo object inside a Data object, based on a geometry
        file.

        Usage:
            data.add_geometry('some_geo_file')
        """
        geom = read_geometry(filename_geo)
        self.geo = geom
        return self

    def filter(self,criterium='latitude',vmin=0,vmax=90):
        """ Method to generate a new data object generated as a subset of data
        from the parent object, filtered along a certain criterium

        Arguments:
            criterium: a string indicating which type of filter is to be applied
                accepted types are : latitude, sza, emission, azimuth,
                local_time, phase, orbit
            vmin: minimal value of the selected criterium
            vmax: maximal value of the selected criterium

        Example :
            NewData = Olddata.filter(criterium='latitude',vmin=0,vmax=45)
        """

        if criterium=='latitude':  # applies filters
            ftr = (self.geo.Lat<vmax) * (self.geo.Lat>vmin)
        elif criterium=='sza':
            ftr = (self.geo.sza<vmax) * (self.geo.sza>vmin)
        elif criterium=='emission':
            ftr = (self.geo.emission<vmax) * (self.geo.emission>vmin)
        elif criterium=='azimuth':
            ftr = (self.geo.azimuth<vmax) * (self.geo.azimuth>vmin)
        elif criterium=='local_time':
            ftr = (self.geo.local_time<vmax) * (self.geo.local_time>vmin)
        elif criterium=='phase':
            ftr = (self.phase<vmax) * (self.phase>vmin)
        elif criterium=='orbit':
            ftr = (self.num_orbit[:,0]>vmin) * (self.num_orbit[:,0]<vmax)
        else:
            raise KeyError('Wrong criterium! Authorised criteria are latitude, '
                           'sza, emission, azimuth, local_time, phase, orbit')

        #creates new filtered objects
        filt_geom = Geo(date = self.geo.date[:,ftr],
                        time = self.geo.time[:,ftr],
                        time_str = self.geo.time_str[ftr,:],
                        L_sun = self.geo.L_sun[ftr],
                        Lat_sc = self.geo.Lat_sc[ftr],
                        Lon_sc = self.geo.Lon_sc[ftr],
                        Lat = self.geo.Lat[ftr],
                        Lon = self.geo.Lon[ftr],
                        Alt_sc = self.geo.Alt_sc[ftr],
                        H = self.geo.H[ftr],
                        phase = self.geo.phase[ftr],
                        sza = self.geo.sza[ftr],
                        emission = self.geo.emission[ftr],
                        local_time = self.geo.local_time[ftr],
                        S_Mdist = self.geo.S_Mdist[ftr],
                        target_dist = self.geo.target_dist[ftr],
                        alt_topo = self.geo.alt_topo[ftr])
        # if radiances measured for several orbits mismatch in shape,
        # the filter might have an issue. This IF is a workaround.
        if self.w0.shape[0] != self.phase.shape[0]:
            filtered = Data(pol = self.pol[:,ftr], dpol = self.dpol[:,ftr],
                            intensity = self.intensity[:,ftr],
                            dI = self.dI[:,ftr],
                            phase = self.geo.phase[ftr], geo = filt_geom,
                            num_orbit = self.num_orbit[ftr,:], w0 = self.w0,
                            w1 = self.w1, wd0 = self.wd0, wd1 = self.wd1,
                            r0 = self.r0, r1 = self.r1,
                            rd0 = self.rd0, rd1 = self.rd1)
        else:
            filtered = Data(pol = self.pol[:,ftr], dpol = self.dpol[:,ftr],
                            intensity = self.intensity[:,ftr],
                            dI = self.dI[:,ftr],
                            phase = self.geo.phase[ftr], geo = filt_geom,
                            num_orbit = self.num_orbit[ftr,0:2],
                            w0 = self.w0[:,ftr], w1 = self.w1[:,ftr],
                            r0 = self.r0[:,ftr], r1 = self.r1[:,ftr],
                            wd0 = self.wd0[:,ftr], wd1 = self.wd1[:,ftr],
                            rd0 = self.rd0[:,ftr], rd1 = self.rd1[:,ftr])



        return filtered


    def split(self):
        """Method to generate a new data object generated with data
        from the parent object but splitted by phase angle.
        New objects are the two branches of phase angle.
        To use when observations have several phase angles in common.

        Usage:
            NewData1, newdata2 = Olddata.split()
        """

        # locate minimum of phase angle
        pos_min_phase = np.argmin(self.phase)

        #creates new filtered object
        ftr = np.arange(0,pos_min_phase) #before minimum
        filt_geom_1 = Geo(date = self.geo.date[:,ftr],
                          time = self.geo.time[:,ftr],
                          time_str = self.geo.time_str[ftr,:],
                          L_sun = self.geo.L_sun[ftr],
                          Lat_sc = self.geo.Lat_sc[ftr],
                          Lon_sc = self.geo.Lon_sc[ftr],
                          Lat = self.geo.Lat[ftr],
                          Lon = self.geo.Lon[ftr],
                          Alt_sc = self.geo.Alt_sc[ftr],
                          H = self.geo.H[ftr],
                          phase = self.geo.phase[ftr],
                          sza = self.geo.sza[ftr],
                          emission = self.geo.emission[ftr],
                          local_time = self.geo.local_time[ftr],
                          S_Mdist = self.geo.S_Mdist[ftr],
                          target_dist = self.geo.target_dist[ftr],
                          alt_topo = self.geo.alt_topo[ftr])
        filtered_1 = Data(pol = self.pol[:,ftr], dpol = self.dpol[:,ftr],
                          intensity = self.intensity[:,ftr],
                          dI = self.dI[:,ftr],
                          phase = self.geo.phase[ftr], geo = filt_geom_1,
                          num_orbit = self.num_orbit[ftr,:],
                          w0 = self.w0[:,ftr], w1 = self.w1[:,ftr],
                          r0 = self.r0[:,ftr], r1 = self.r1[:,ftr],
                          wd0 = self.wd0[:,ftr], wd1 = self.wd1[:,ftr],
                          rd0 = self.rd0[:,ftr], rd1 = self.rd1[:,ftr])

        ftr = np.arange(pos_min_phase,len(self.phase)) #after minimum
        filt_geom_2 = Geo(date = self.geo.date[:,ftr],
                          time = self.geo.time[:,ftr],
                          time_str = self.geo.time_str[ftr,:],
                          L_sun = self.geo.L_sun[ftr],
                          Lat_sc = self.geo.Lat_sc[ftr],
                          Lon_sc = self.geo.Lon_sc[ftr],
                          Lat = self.geo.Lat[ftr],
                          Lon = self.geo.Lon[ftr],
                          Alt_sc = self.geo.Alt_sc[ftr],
                          H = self.geo.H[ftr],
                          phase = self.geo.phase[ftr],
                          sza = self.geo.sza[ftr],
                          emission = self.geo.emission[ftr],
                          local_time = self.geo.local_time[ftr],
                          S_Mdist = self.geo.S_Mdist[ftr],
                          target_dist = self.geo.target_dist[ftr],
                          alt_topo = self.geo.alt_topo[ftr])
        filtered_2 = Data(pol=self.pol[:,ftr], dpol=self.dpol[:,ftr],
                          intensity=self.intensity[:,ftr],
                          dI = self.dI[:,ftr],
                          phase=self.geo.phase[ftr], geo=filt_geom_2,
                          num_orbit=self.num_orbit[ftr,0:2],
                          w0 = self.w0[:,ftr], w1 = self.w1[:,ftr],
                          r0 = self.r0[:,ftr], r1 = self.r1[:,ftr],
                          wd0 = self.wd0[:,ftr], wd1 = self.wd1[:,ftr],
                          rd0 = self.rd0[:,ftr], rd1 = self.rd1[:,ftr])


        return filtered_1, filtered_2


    def mean_data_phase(self,step=1,wvl=1101,band=False):
        """Reduction of data's noise in phase angle at a given wavelength

        Arguments:
            step: step
            wvl: wavelength on which the average must be done
            band: if True, the average is done in the spectral band. If False,
                only DOTS are considered.
        Output:
            adds new attributes to the object:
                x_int: abscissae of the averaged points
                pol_mean: averaged polarization
                Spol_mean: error bar of the averaged polarization


        xbins : width of meantime
        n : number of points
        pol_mean : data averaged in phase
        Spol_mean : errorbar of data averaged
        w_i : statistical weight"""

        # Initialisation
        if band==False:  #Continuum
            phase = self.phase
            sza = self.geo.sza
            emi = self.geo.emission
            pol = self.pol
            dpol = self.dpol

            V = self.wd0[:,0]

        elif band==True:  #In the band
            self = self.get_bande()

            phase = self.phase
            pol = self.band_pol
            dpol = self.band_dpol

            V = self.band_wd0[:,0]

        # Check if user wvl is available
        # in particular if looking in the band
        if (wvl<min(V)) or (wvl>max(V)):
            print("Warning ! wvl value is out of range")

        amin = min(phase)
        amax = max(phase)
        xbins = np.arange(amin,amax,step)

        #Searching for the nearest spectel from wvl
        wdx = np.where(abs(V-wvl) == np.nanmin(abs(V-wvl)))[0]
        print("Retrieved wavelength : {:4f} nm".format(V[wdx][0]))

        #Creation of final arrays
        pol_mean = np.zeros(len(xbins)-1)
        sza_mean = np.zeros(len(xbins)-1)
        emi_mean = np.zeros(len(xbins)-1)
        Spol_mean = np.zeros(len(xbins)-1)

        #Computation of abscissa
        x_int = xbins[:-1] + np.diff(xbins)/2.

        # Loop on bins
        for i in np.arange(len(xbins)-1):
            x_indices = np.where( (phase>xbins[i]) & (phase<xbins[i+1]) )
                #mask of indices for each meantime
            x_indices = x_indices[0]
            n = len(x_indices)

            #Computation of points averaged in phase
            w_i = 1/(dpol[wdx,x_indices])**2
            w_sum = sum(w_i)  #sum of statistical weights

            pol_i = w_i*pol[wdx,x_indices]
            pol_sum = sum(pol_i)  #sum of weightened polarisations
            # locating the phase angle closest to bin phase
            phase_eq = np.where(abs(phase-x_int[i]) == np.nanmin(abs(phase-x_int[i])))[0]
            print(phase_eq)
            phase_eq = phase_eq[:1]
            pol_mean[i] = pol_sum/w_sum
            sza_mean[i] = sza[phase_eq]
            emi_mean[i] = emi[phase_eq]

            #Computation of errorbars
            Spol_i = w_i*((pol[wdx,x_indices]-pol_mean[i])**2)
            Spol_sum = sum(Spol_i)
            Spol_mean[i] = np.sqrt(Spol_sum/(n-1)/w_sum)

        self.x_int = x_int
        self.sza_mean = sza_mean
        self.emi_mean = emi_mean
        self.pol_mean = pol_mean
        self.Spol_mean = Spol_mean

        return self

    def get_bande(self,lw_noise='lw_NEB_g2_t28.txt',sw_noise='sw_NEB_g2_t28.txt'):
        """This method returns the data for all wavelength in the CO2 band

        Arguments:
            lw_noise: NEB of detectors in the LW range
            sw_noise: NEB of detectors in the SW range


        """

        #Initialisation
        wl_d0 = self.w0
        wl_d1 = self.w1
        radf_d0 = self.r0
        radf_d1 = self.r1

        nobs = len(wl_d0[0,:])

        r0_invalid =  (radf_d0<0) # identification of negative radiance values
        r1_invalid =  (radf_d1<0)

        radf_d0[r0_invalid] = np.nan # negative values are set to NaN
        radf_d1[r1_invalid] = np.nan

        #Noise equivalent brightness
        neb_sw = np.genfromtxt(sw_noise,skip_header=1)
        neb_lw = np.genfromtxt(lw_noise,skip_header=1)

        neb = np.vstack((neb_lw,neb_sw))

        #def of band start and end in nanometers
        start_band = 1246.
        end_band = 1490.

        # identifying number of spectels
        indices0 = np.where( (wl_d0[:,0]>start_band) & (wl_d0[:,0]<end_band)  )[0]
        nspectels = len(indices0)

        # creation of final arrays
        band_wvl0 = np.zeros([nspectels,nobs]) # shape(nb of wvl, number of measures)
        band_wvl1 = np.zeros([nspectels,nobs])

        band_rad0 = np.zeros([nspectels,nobs])
        band_rad1 = np.zeros([nspectels,nobs])

        neb_det0 = np.zeros([nspectels])
        neb_det1 = np.zeros([nspectels])
        dpol = np.zeros([nspectels,nobs])

        #selection of wavelengths in the NEB data
        noise_idx = (neb[:,0]>start_band) * (neb[:,0]<end_band)
        neb_wvl = neb[noise_idx,0]

        #mask on the NEB data  so we get the noise on detectors for all
        # reference wavelengths
        neb_det0 = neb[noise_idx,1]
        neb_det1 = neb[noise_idx,2]

        filtre = np.argsort(neb_wvl)
        neb_wvl = neb_wvl[filtre]
        neb_det0 = neb_det0[filtre]
        neb_det1 = neb_det1[filtre]


        # Getting the data for all wavelength in the band and for all
        # observation points
        for j in range(nobs):  #loop on observation times

            # identification of wavelengths by
            # using a mask (excluding nans on radiance))
            indices0 = np.where( (wl_d0[:,j]>start_band)
                                & (wl_d0[:,j]<end_band) )[0]
            indices1 = np.where( (wl_d1[:,j]>start_band)
                                & (wl_d1[:,j]<end_band) )[0]

            if (len(indices0)==0) or  (len(indices1)==0):
                band_wvl0[:,j] = np.nan
                band_wvl1[:,j] = np.nan

                band_rad0[:,j] = np.nan
                band_rad1[:,j] = np.nan
            else:
                band_wvl0[:,j] = wl_d0[indices0,j]
                band_wvl1[:,j] = wl_d1[indices1,j]

                band_rad0[:,j] = radf_d0[indices0,j]
                band_rad1[:,j] = radf_d1[indices1,j]

        pol = (band_rad1 - band_rad0) / (band_rad0 + band_rad1)

        neb_0 = np.interp(band_wvl0, neb_wvl, neb_det0)
        neb_1 = np.interp(band_wvl1, neb_wvl, neb_det1)

        #calculus of the error on polarization for all observation points
        for i,d in enumerate(band_rad0[0,:]):
            dpol[:,i] = np.sqrt( ((2*band_rad0[:,i]*neb_1[:,i])
                                  /(band_rad0[:,i] + band_rad1[:,i])**2)**2
                                + ( (2*band_rad1[:,i]*neb_0[:,i])/
                                   (band_rad0[:,i] + band_rad1[:,i])**2)**2 )

        intensity = (band_rad0 + band_rad1)

        self.band_wd0 = band_wvl0
        self.band_wd1 = band_wvl1
        self.band_rd0 = band_rad0
        self.band_rd1 = band_rad1
        self.band_pol = pol
        self.band_dpol = dpol
        self.band_intensity = intensity

        return self

    def show_band(self, vmin=-10, vmax=10):
        """ Displays the data cube related to the absorption band"""
        fig = mpl.figure()
        ax = fig.add_subplot(111)
        cube = ax.pcolor(self.phase, self.band_wd0[:,0],
                         100*self.band_pol,
                         vmin=vmin, vmax=vmax, cmap=mpl.cm.YlOrRd)
        ax.set_xlabel('Phase angle', size=14)
        ax.set_ylabel('Wavelength (nm)', size=14)
        fig.tight_layout()
        cb = mpl.colorbar(cube, pad=0.02)
        cb.set_label('Degree of linear polarization', size=14)

    def phase_pol_colat(self,wvl_idx,font_size=14):
        wvl_list = [0.650,0.760,0.855,0.980,1.101,1.160,1.198,1.242,1.274,1.305,1.324,1.515,1.553,1.625]
        numorb0 = int(self.num_orbit[0,0])  # orbit numbers for legend
        numorb1 = int(self.num_orbit[0,1])
        wvl = wvl_list[wvl_idx]
        fig = mpl.figure(figsize=(20,10))
        for window in np.arange(9):
            wvl_idx = window+4
            wvl = wvl_list[wvl_idx]
            ax = fig.add_subplot(3,4,window+1)
            l1 = ax.errorbar(self.geo.phase,100*self.pol[wvl_idx,:],yerr=100*self.dpol[wvl_idx,:],fmt='.', label='Orbit {:04d}-{:02d} @{:04d}nm'.format(numorb0,numorb1,int(wvl*1000)))  # plot data
            ax2 = ax.twiny()  # create twin axis
            l2 = ax2.plot(self.geo.phase,100*self.pol[wvl_idx,:],alpha=0)  # plot data but invisible!
            xticks = ax.get_xticks()  # get ticks of bottom axis
            xticks = xticks[1::2]
            nticks = len(xticks)
            xticks2 = np.zeros(nticks)
            for i,t in enumerate(xticks):  # for each tick
                idx = np.argmin(abs(self.geo.phase-t))  # find index that is closest
                xticks2[i] = round(self.geo.Lat[idx],2)  # get value of latitude

            ax2.set_xticks(xticks)  # set top ticks to same position
            ax2.set_xticklabels(xticks2)  # give them values
            #ax2.set_xlabel('latitude',size=8)
            ax.set_xlabel('Phase angle',size=8)
            #ax.legend(loc='lower right')
            #ax.set_ylabel('Degree of linear polarization (%)',size=font_size)
        fig.suptitle('Orbit {:04d}-{:02d} @{:04d}nm'.format(numorb0,numorb1,int(wvl*1000)))
        fig.tight_layout(pad=1.5)
        mpl.subplots_adjust(top=0.5)
        mpl.savefig('phase_pol_colat_{:04d}.png'.format(numorb0,int(1000*wvl)))


#----------------------
# FUNCTIONS DEFINITION
#----------------------

def read_radiance(filename):
    """ READ_RADIANCE
    This function reads the 0B level data from SPICAV instrument on VEx.

    INPUT :
    --- filename : name of the data file to be read.

    OUTPUTS :
        data: a Data object with basic attributes.

    Note:
        wl_d0 : wavelength array for detector 0. Has shape (nb of wavelenghts,
            nb of observations).
        wl_d1 : idem, but for detector 1.
        radf_d0 : radiance array for dectector 0. Has shape (nb of wavelenghts,
            nb of observations).
        radf_d1: radiance array for dectector 1. Has shape (nb of wavelenghts,
            nb of observations).
        time_data : ?

    EXAMPLE:
        data = read_radiance('SPIV_0BR_0123A02_N_polar.DAT')

    COMPLEMENTS :
    # Variables :
    #  - npt : number of spectral points
    #  - fnb : number of spectra
    #  - channels : number of channels
    #  - nbl : number of blocks
    #  - ngt : index, if ngt=1, night observation, if ngt=0, day obs.
    #
    #  - wl_d0 : wavelength vector for channel 0
    #  - radf_d0 : radiance for channel 0 in W/m2/um/sr
    #  - tok_d0 : dark current calibration for channel 0, already taken into
    #  account in radiance

    """

    try:
        fid = open(filename,'rb')
    except IOError:
        print('Error, file {} does not exist. '.format(filename))

    # reads parameters for the next steps
    s = np.fromfile(fid,count=5,dtype='int16')
    channels = s[0] #number of detectors
    npt = s[1] # number of spectral points
    fnb = s[2] # number of spectra
    nbl = s[3] # number of blocks
    ngt = s[4] # index, if ngt=1 it is night obs. if ngt=0, it is day obs.

    #creation of tables
    wl_d0 = np.zeros([npt,fnb])
    wl_d1 = np.zeros([npt,fnb])
    tok_d0 = np.zeros([npt,fnb])
    tok_d1 = np.zeros([npt,fnb])
    clb_d0 = np.zeros([npt,fnb])
    clb_d1 = np.zeros([npt,fnb])
    radf_d0 = np.zeros([npt,fnb])
    radf_d1 = np.zeros([npt,fnb])

    #Reading the data
    for i in range(0,fnb-1):
        time_dat = np.fromfile(fid,count=6,dtype='int16')
        tmp = np.fromfile(fid,count=1,dtype='double')
        time_dat = np.append(time_dat,tmp)
        wl_d0[:,i] = np.fromfile(fid,count=npt,dtype='float32')
        wl_d1[:,i] = np.fromfile(fid,count=npt,dtype='float32')
        tok_d0[:,i] = np.fromfile(fid,count=npt,dtype='float32')
        tok_d1[:,i] = np.fromfile(fid,count=npt,dtype='float32')
        clb_d0[:,i] = np.fromfile(fid,count=npt,dtype='float32')
        clb_d1[:,i] = np.fromfile(fid,count=npt,dtype='float32')
        radf_d0[:,i] = np.fromfile(fid,count=npt,dtype='float32')
        radf_d1[:,i] = np.fromfile(fid,count=npt,dtype='float32')

        if (i==0):
            time_data = time_dat

    fid.close()

    print("nbl =", nbl)
    print("fnb =", fnb)

    data = Data()
    data.w0 = wl_d0
    data.w1 = wl_d1
    data.r0 = radf_d0
    data.r1 = radf_d1
    data.time = time_data
    data.filename = filename

    return data

def read_geometry(filename):

    """READ_GEOMETRY
    This function reads the geometry file of name *filename*

    INPUT :
    --- filename : name of the file containing the geometric parameters.
    OUTPUT :
    --- a Geometric object containing all geometric parameters
    see the Geo class for details.
    """

    if os.path.isfile(filename)==0:
        print("Geometry file {} does not exist.".format(filename))

    header_length = 1
    geodata = np.genfromtxt(filename,dtype='str',skip_header=header_length)
    #reads text file, skipping the header

    date = np.vstack((geodata[:,0],geodata[:,1],geodata[:,2])) #reading date
    time0 = geodata[:,3]  #reading the time string
    time_str = []
    for i,zz in enumerate(date[0,:]):
        time_str.append(' '.join(date[:,i]) +' ' + time0[i])
        #makes a date + time string
    time_str = np.asarray(time_str, dtype='c')
    time_geo_hour = np.zeros(0)
    time_geo_minute = np.zeros(0)
    time_geo_second = np.zeros(0)
    for temp in time0:
        time_geo_hour = np.append(time_geo_hour,
                                  float(temp[0:2]))  #reading hours in it
        time_geo_minute = np.append(time_geo_minute,
                                    float(temp[3:5]))  #reading minutes in it
        time_geo_second = np.append(time_geo_second,
                                    float(temp[6:11])) # reading seconds in it

    # Storing data
    time = np.vstack((time_geo_hour,time_geo_minute,time_geo_second))
    L_sun = (geodata[:,4]).astype(float)
    Lat_sc = (geodata[:,5]).astype(float)
    Lon_sc = (geodata[:,6]).astype(float)
    Lat = (geodata[:,7]).astype(float)
    Lon = (geodata[:,8]).astype(float)
    Alt_sc = (geodata[:,9]).astype(float)
    H = (geodata[:,10]).astype(float)
    PHASE = (geodata[:,11]).astype(float)
    SZA = (geodata[:,12]).astype(float)
    EMISSN = (geodata[:,13]).astype(float)

    ltime0 = geodata[:,14]  #time string
    ltime = np.zeros(0)
    for temp in ltime0:
        ltime = np.append(ltime,
                          float(temp[0:2])
                          + float(temp[3:5])/60.
                          + float(temp[6:11])/3600.) # getting minutes out of it

    S_Mdist = (geodata[:,15]).astype(float)
    TARGETdist = (geodata[:,16]).astype(float)
    Alt_topo = (geodata[:,17]).astype(float)

    #Creation of geometric object
    geom = Geo(date = date,
               time = time,
               time_str = time_str,
               L_sun = L_sun,
               Lat_sc = Lat_sc,
               Lon_sc = Lon_sc, Lat = Lat,
               Lon = Lon, Alt_sc = Alt_sc,
               H = H, phase = PHASE,
               sza = SZA, emission = EMISSN,
               local_time = ltime, S_Mdist = S_Mdist,
               target_dist = TARGETdist, alt_topo = Alt_topo)

    return geom

def clean_nan(in_array):
    """Cleans an array to make it contain only finite (not NaN) values)"""
    not_nan = np.isfinite(in_array)
    out_array = in_array[not_nan]
    return out_array


def processing_orbits(orbit_n,orbit_a,orbit_max=200,list=False,
                      orbit_range=False,
                      path_data='/Users/gouravmahapatr/Dropbox/PhD/spicav_data/calibrated_data2013/',
                      calibration_year=2013,
                      path_geo='/Users/gouravmahapatr/Dropbox/PhD/spicav_data/geom/',
                      lw_noise='/Users/gouravmahapatr/Dropbox/PhD/spicav_data/lw_NEB_g2_t28.txt',
                      sw_noise='/Users/gouravmahapatr/Dropbox/PhD/spicav_data/sw_NEB_g2_t28.txt'):
    """
    Processing_orbits
    Processes orbits and gives useful outputs

    This is a convienience routine. For a given (set of) orbit(s), it
    reads the data, the geometry, treats the radiance for the
    polarization measure, and stores all the results in a Data object.

    COMPULSORY INPUTS :
        orbit_n : Orbit number
        orbit_a : orbit sequence number
    OPTIONNAL INPUTS :
        orbit_range : if this parameter is set to true, the parameter orbit_max
            is used, and all available orbits from orbit_n to orbit_max are
            combined.
        orbit_max : maximal orbit if orbit_range is set to True
        list : if set to True, the function can take vectors as orbit_n and
        orbit_a. All indicated orbits are combined.
        path_data : path of the data files.

    OUTPUT:
        a Data object
    """

    index = 'N'
    output = {}
    counter = 0

    # SINGLE ORBIT
    if list==False and orbit_range==False:

        filename = (path_data
                    +'SPIV_0BR_'
                    + str('{:04d}'.format(orbit_n))+ 'A'
                    + str('{:02d}'.format(orbit_a))+ '_'
                    + index + '_05_calibr_'+str(calibration_year)
                    +'_polar_virtis.DAT')
        filename_geo = (path_geo
                        +'SPIV_IR_' + str('{:04d}'.format(orbit_n))
                        + 'A' + str('{:02d}'.format(orbit_a))
                        + '_' + index + '_2012.txt')

        if os.path.isfile(filename)==False:
            raise IOError('Error'+filename+' is not an existing file')
        else:
            print(filename)
            data = Data()

            data = read_radiance(filename)  #read data
            data.adjust_radiance(lw_noise=lw_noise,sw_noise=sw_noise)  #average

            # Reading geometry
            data.add_geometry(filename_geo)
            lenn = len(data.geo.phase)
            num_orbit = np.zeros((lenn,2))
            num_orbit[:,0:2] = np.array([orbit_n,orbit_a])

            data.phase = data.geo.phase
            data.num_orbit = num_orbit
            data_final = data
            counter += 1

    # MULTIPLE ORBITS
    elif list==True:

        #check of vector shape
        if len(np.shape(orbit_n)) != 1 | len(np.shape(orbit_a)) !=1:
            print("Error, wrong vector shape !")

        counter=0

        data_final = Data()

        for i,num in enumerate(orbit_n):  #read all mentioned orbits
            ind = orbit_a[i]
            filename = (path_data+'SPIV_0BR_'
                        + str('{:04d}'.format(num))+ 'A'
                        + str('{:02d}'.format(ind))+ '_'
                        + index + '_05_calibr_'
                        +str(calibration_year)+'_polar_virtis.DAT')
            filename_geo = (path_geo+'SPIV_IR_'
                            + str('{:04d}'.format(num))
                            + 'A' + str('{:02d}'.format(ind))
                            + '_' + index + '_2012.txt')

            if os.path.isfile(filename)==False :
                print('absent file :'+filename)
            elif os.path.isfile(filename_geo)==False:
                print('absent file :'+filename_geo)
            else:
                data = read_radiance(filename)  #read data
                data.adjust_radiance(lw_noise=lw_noise,
                                     sw_noise=sw_noise)  #average
                data.add_geometry(filename_geo)
                lenn = len(data.geo.phase)
                num_orbit = np.zeros((lenn,2))
                num_orbit[:,0:2] = np.array([num,ind])

                data.phase = data.geo.phase
                data.num_orbit = num_orbit

                data_final += data
                counter += 1


    # ORBIT RANGE
    elif orbit_range==True:
        counter=0
        data_final = Data()

        for num in range(orbit_n,orbit_max):
            for ind in range(0,20):
                filename = (path_data
                            +'SPIV_0BR_' + str('{:04d}'.format(num))
                            + 'A' + str('{:02d}'.format(ind))
                            + '_' + index + '_05_calibr_'
                            + str(calibration_year)+'_polar_virtis.DAT')
                filename_geo = (path_geo
                                +'SPIV_IR_'
                                + str('{:04d}'.format(num))+ 'A'
                                + str('{:02d}'.format(ind))+ '_'
                                + index + '_2012.txt')

                if os.path.isfile(filename)==False :
                    print('absent file :'+filename)
                elif os.path.isfile(filename_geo)==False:
                    print('absent file :'+filename_geo)
                else:
                    data = read_radiance(filename)  #read data
                    data.adjust_radiance(lw_noise=lw_noise,
                                         sw_noise=sw_noise)  #average
                    data.add_geometry(filename_geo)
                    lenn = len(data.geo.phase)
                    num_orbit = np.zeros((lenn,2))
                    num_orbit[:,0:2] = np.array([num,ind])

                    data.phase = data.geo.phase
                    data.num_orbit = num_orbit

                    data_final += data
                    counter += 1

    print('{} file(s) have been read'.format(counter) )

    return data_final


def average_curve(data,idx,n_bins=90,disk=0,step=1,med=1):
    """computes average curve of polarization as a function of phase angle
    Call : x,y = average_curve(data,wvl number)
    options :
        disk: disk integration
        step: step between bin steps
        med: if 1 returns median instead of average. If 0 returns average
    """
    x = data.phase
    y = data.pol[idx,:]
    rad = data.intensity[idx,:]

    gval = np.isfinite(y) #good values not NaN
    y = y[gval]
    x = x[gval]
    rad = rad[gval]
    bins = np.arange(0,90,step)  # set custom bins
    medians = np.zeros(len(bins)-1) # set a vector of medians

    for i in np.arange(len(bins)-1):
        medians[i] = np.median(y[(x>bins[i])*(x<bins[i+1])]) # compute media of values inside bin

    if disk==0:
        n,av_x = np.histogram(x,bins=bins)
        av_x = av_x[:-1] + 0.5*np.diff(av_x)
        sum_y = np.histogram(x,bins=bins,weights=y)[0]
        av_pol = sum_y/n
    elif disk==1:
        n,av_x = np.histogram(x,bins=bins)
        av_x = av_x[:-1] + 0.5*np.diff(av_x)
        sum_y = np.histogram(x,bins=bins,weights=y*rad)[0]
        sum_r = np.histogram(x,bins=bins,weights=rad)[0]
        av_pol = sum_y/sum_r

    if med==1:
        return  av_x,medians
    elif med==0:
        return av_x,av_pol


def temporal_bins(data,vmin=400,vmax=2300,step=100):

    bins = np.arange(vmin,vmax,step)
    nbins = len(bins)
    for i in np.arange(nbins-1):
        band = data.filter(criterium='orbit',vmin=bins[i],vmax=bins[i+1])
        x,y = average_curve(band,4,step=2)
        mpl.plot(x,100*y,label='{}-{}'.format(bins[i],bins[i+1]))
        #mpl.hexbin(band.phase,100*band.pol[4,:],mincnt=5)
        mpl.xlim(0,90)
        mpl.ylim(-8,6)
        mpl.savefig('band_{}'.format(i))
        #mpl.clf()

def show_path(data,ax=None,nstep=10,png=False,LT=False, legendpos='upper left', cvar=None):
    """ This function allows to see the path of the spacecraft on the planet
    and the FOV
    INPUT:
        - data : a Data object
        - nstep : (optionnal) the step between different points plotted
    OUTPUT:
        a graphical output
    """
    # Defining vectors of data to plot
    # "sc" meaning spacecraft
    xsc = data.geo.Lon_sc[::nstep]
    ysc = data.geo.Lat_sc[::nstep]
    x = data.geo.Lon[::nstep]
    y = data.geo.Lat[::nstep]
    # Getting position and illumation angles on the edge of FOV
    lats, longs, phase, sza, emi, beta, dist, local_time = fovint(data.geo.time_str)
    # If the FOV does not intersect the planet, put a NaN
    emi[lats==-1000] = np.nan
    lats[lats==-1000] = np.nan
    longs[longs==-1000] = np.nan
    phase[phase==-1000] = np.nan
    local_time[local_time==-1000] = np.nan

    # Boundaries of plotable values
    firstlong = np.where(~np.isnan(np.sum(longs,axis=0)))[0][0]
    firstlat = np.where(~np.isnan(np.sum(lats,axis=0)))[0][0]
    firstlt = np.where(~np.isnan(np.sum(local_time,axis=0)))[0][0]
    lastlong = np.where(~np.isnan(np.sum(longs,axis=0)))[0][-1]
    lastlat = np.where(~np.isnan(np.sum(lats,axis=0)))[0][-1]
    lastlt = np.where(~np.isnan(np.sum(local_time,axis=0)))[0][-1]
    firstphase = np.max([firstlong,firstlat,firstlt])
    lastphase = np.max([lastlong,lastlat,lastlt])
    maxx = np.nanmax(longs[:,0])
    miny = np.nanmin(lats[:,0])
    maxy = np.nanmax(lats[:,0])
    minphase = np.min(phase[0,:])
    idxminphase = np.where(phase[0,:]==minphase)[0]
    #leftx = np.where(longs[:,0]==minx)[0]
    #rightx = np.where(longs[:,0]==maxx)[0]
    #leftx = np.where(longs[:,0]==minx)[0]
    #rightx = np.where(longs[:,0]==maxx)[0]

    if ax == None:
        fig = mpl.figure()
        ax = fig.add_subplot(111)

    if LT == False:
        minlongs = map(np.nanmin,longs.T)
        maxlongs = map(np.nanmax,longs.T)
        minlats = map(np.nanmin,lats.T)
        maxlats = map(np.nanmax,lats.T)
        Xs = np.append(minlongs,maxlongs[::-1])
        Ys = np.append(minlats,maxlats[::-1])
        #ax.fill(Xs,Ys,fill=False)
        ax.fill(longs[1:,firstlong],lats[1:,firstlat],lw=1, label='FOV at start',fill=False)
        #ax.scatter(longs[0,0], lats[0,0], lw=0, c='g', label='center at start')
        ax.fill(longs[1:,lastlong],lats[1:,lastlat],fill=False, lw=2, label='FOV at end')
        #ax.scatter(longs[0,-1], lats[0,-1], lw=0, c='r',label='center at end')
        if cvar=='emi':
            sc = ax.scatter(longs[0,firstphase:lastphase],lats[0,firstphase:lastphase],c=emi[0,firstphase:lastphase],lw=0,s=100, cmap='copper')
        else:
            sc = ax.scatter(longs[0,firstphase:lastphase],lats[0,firstphase:lastphase],c=phase[0,firstphase:lastphase],lw=0,s=100, cmap='copper')
        cb = fig.colorbar(sc,pad=0.02)
        ax.set_xlabel('Longitude', size=14)

    elif LT == True:
        minlongs = map(np.nanmin,local_time.T)
        maxlongs = map(np.nanmax,local_time.T)
        minlats = map(np.nanmin,lats.T)
        maxlats = map(np.nanmax,lats.T)
        Xs = np.append(minlongs,maxlongs[:])
        Ys = np.append(minlats,maxlats[:])
        #ax.fill(Xs,Ys,fill=False)
        ax.fill(local_time[1:,firstlt],lats[1:,firstlat],lw=1, label='FOV at start',fill=False)
        #ax.scatter(local_time[0,0], lats[0,0], lw=0, c='g', label='center at start')
        #ax.fill(local_time[1:,idxminphase],lats[1:,idxminphase],fill=False, lw=2, label='FOV at minimum phase')
        ax.fill(local_time[1:,lastlt],lats[1:,lastlat],fill=False, lw=2, label='FOV at end')
        #ax.scatter(local_time[0,-1], lats[0,-1], lw=0, c='r',label='center at end')
        if cvar=='emi':
            sc = ax.scatter(local_time[0,firstphase:lastphase],lats[0,firstphase:lastphase],c=emi[0,firstphase:lastphase],lw=0,s=100, cmap='copper')
        else:
            sc = ax.scatter(local_time[0,firstphase:lastphase],lats[0,firstphase:lastphase],c=phase[0,firstphase:lastphase],lw=0,s=100, cmap='copper')
        cb = fig.colorbar(sc,pad=0.02)
        ax.set_xlabel('Local time', size=14)


    if cvar=='emi':
        cb.set_label('Emission angle')
    else:
        cb.set_label('Phase angle')
    ax.legend(loc=legendpos)
    ax.set_ylabel('Latitude (N)', size=14)
    ax.set_title('FOV and trajectory of SPICAV Orbit {:4.0f}-{:2.0f}'.format(data.num_orbit[0,0],data.num_orbit[0,1]))
    fig.tight_layout()
    #mpl.savefig('path_fov_tmp_{:4d}-{:1d}.png'.format(data.num_orbit[0,0].astype('int'),data.num_orbit[0,1].astype('int')))
    mpl.savefig('path_fov_tmp_{:3d}-{:1d}.pdf'.format(data.num_orbit[0,0].astype('int'),data.num_orbit[0,1].astype('int')))

    return ax, Xs, Ys

def calc_dist(data,R):
    """ computes the distance between two points given by their lat and
    longitude on a sphere with radius R."""

    lats, longs, phase, sza, emi, beta, dist, local_time = fovint(data.geo.time_str)
    lats[lats==-1000] = np.nan
    longs[longs==-1000] = np.nan
    phase[phase==-1000] = np.nan
    local_time[local_time==-1000] = np.nan

    nans = np.isnan(lats)[0]
    okaylats = np.where(~np.isnan(lats))[1]
    okaylongs = np.where(~np.isnan(longs))[1]

    latA = np.radians(okaylats[0])
    latB = np.radians(okaylats[-1])
    longA = np.radians(okaylongs[0])
    longB = np.radians(okaylongs[-1])

    dl = longB - longA
    S = np.arccos(np.sin(latA)*np.sin(latB) + np.cos(latA)*np.cos(latB)*np.cos(dl))
    dist = S*R

    return dist

def find_orbits(omin=None,omax=None,orb_type=None):
    ''' This function scans through all the available orbit files
    from the geometry folder of SPICAV and obtains a list of orbits,
    with or without given constrains.
    
    Optional inputs:
        omin: Minimum orbit number to be read
        omax: Maximum orbit number to be read 
    
    orb_type:
        'None': All the orbits that exist are given in an array. 
        
        
    Author: G. Mahapatra
        '''
        
    orbit_list = []
    
    path_data = '/Users/gouravmahapatr/Dropbox/PhD/spicav_data/calibrated_data2013/'
    path_geo = '/Users/gouravmahapatr/Dropbox/PhD/spicav_data/geom/'
    
    index = 'N' # no 'A' data available yet
    
    calibration_year = 2013
    
    if omin==None and omax==None: 
        omin = 1
        omax = 3200
    
    for i in range(omin,omax):
        orbit_n = i
        for j in range(1,25):
            orbit_a = j
            filename = (path_data
                            +'SPIV_0BR_'
                            + str('{:04d}'.format(orbit_n))+ 'A'
                            + str('{:02d}'.format(orbit_a))+ '_'
                            + index + '_05_calibr_'+str(calibration_year)
                            +'_polar_virtis.DAT')
            filename_geo = (path_geo
                            +'SPIV_IR_' + str('{:04d}'.format(orbit_n))
                            + 'A' + str('{:02d}'.format(orbit_a))
                            + '_' + index + '_2012.txt')
        
            if os.path.isfile(filename)==False:
                #raise IOError('Error'+filename+' is not an existing file')
                continue
            else:
                orbit_list.append((orbit_n,orbit_a))
        
    return np.array(orbit_list)


def getSPIRflux(lmin=None, lmax=None,idx=0,det=0):
    """
    This function gets the NORMALIZED flux vectors from the SPICAV-IR routine
    within a requested given range of wavelengths.
    You can input an orbital epoch also, default is set at index=0.
    
    Returns: radiance at a given detector,wavelength
    
    Authr: G. Mahapatra
    """
    if det == 0:
        r0 = data.r0[:,idx]
        w0 = data.w0[:,idx]
    else:
        r0 = data.r1[:,idx]
        w0 = data.w1[:,idx]
        
    # drop the "nan's" and normalize the radiance
    w0 = w0[~np.isnan(r0)]
    r0 = r0[~np.isnan(r0)]  # remove nan's
        
    if (lmin==None) & (lmax == None):
        # get the SPICAV-IR radiance and wavelength data
        r = r0
        w = w0
    else:
        r = r0[(w0>lmin)&(w0<lmax)]
        w = w0[(w0>lmin)&(w0<lmax)]
        
    r = r/r[-1]
        
    return r,w
    
def getSPIRpol(lmin=None, lmax=None,idx=0,det=0,orbn=1478,orba=8):  
    """
    This function returns polarized flux data for a given orbit geometry =idx.
    
    It does not do any processing as done by Loic Rossi to arrive at his selected 
    bands.
    
    pol = (det1 - det0)/(det0 + det1)
    
    Returns: pol,wavelength
    
    Authr: G. Mahapatra
    
    """
    
    data = processing_orbits(orbn,orba)
    
    r0 = data.r0[:,idx]
    w0 = data.w0[:,idx]
    
    r1 = data.r1[:,idx]
    w1 = data.w1[:,idx]
    
    
    if (lmin==None) & (lmax == None):
    # get the SPICAV-IR radiance and wavelength data
        r0 = r0
        w0 = w0
        r1 = r1
        w1 = w1
    else:
        lidx = (w0>lmin)&(w0<lmax)
        r0 = r0[lidx]
        w0 = w0[lidx]
        r1 = r1[lidx]
        w1 = w1[lidx]
        
    # calculate the polarization at all the available geometries
    # polarization is calculated as (det1 - det0)/(det0 + det1)
    pol = np.nan_to_num((r1 - r0)/(r0 + r1))
    
    return pol,w1
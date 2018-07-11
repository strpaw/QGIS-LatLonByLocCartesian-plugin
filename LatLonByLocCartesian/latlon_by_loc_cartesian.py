# -*- coding: utf-8 -*-
"""
/***************************************************************************
 LatLonByLocCartesian
                                 A QGIS plugin
 latitude, longitude calcucation based on local cartesian coordinate system
                              -------------------
        begin                : 2018-07-09
        git sha              : $Format:%H$
        copyright            : (C) 2018 by Paweł Strzelewicz
        email                : pawel.strzelewicz@wp.pl
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
"""
from PyQt4.QtCore import *
from PyQt4.QtGui import QAction, QIcon, QFileDialog, QMessageBox, QWidget
from qgis.core import *
# Initialize Qt resources from file resources.py
import resources
# Import the code for the dialog
from latlon_by_loc_cartesian_dialog import LatLonByLocCartesianDialog
import os.path

import re
import math
import csv
import datetime

# Constants

# Parameters of WGS84 ellipsoid
WGS84_A = 6378137.0         # semi-major axis of the WGS84 ellipsoid in m
WGS84_B = 6356752.314245    # semi-minor axis of the WGS84 ellipsoid in m
WGS84_F = 1 / 298.257223563 # flatttening of the WGS84 ellipsoid

# Special constants to use instead of False, to avoid ambigous where result of function might equal 0 and 
# and result of fucntion will be used in if statements etc.
VALID = 'VALID'
NOT_VALID = 'NOT_VALID'

# Units of measure
UOM_M = 'M'     # meters
UOM_KM = 'KM'    # kilometers
UOM_NM = 'NM'    # nautical miles
UOM_FEET = 'FEET'  # feet
UOM_SM = 'SM'    # statue miles

# Value types constant
V_AZM = 'AZM'
V_MAG_VAR = 'MAG_VAR'
V_LAT = 'LAT'
V_LON = 'LON'

# DMS, DM format separators, e.g. N32-44-55.21, N32 44 55.21, N32DEG 44 55.21
S_SPACE = ' '     # Blank, space separator
S_HYPHEN = '-'     # Hyphen separator
S_WORD_DEG = 'DEG'
S_WORD_MIN = 'MIN'
S_WORD_SEC = 'SEC'
S_ALL = [S_SPACE, S_HYPHEN, S_WORD_DEG, S_WORD_MIN, S_WORD_SEC]

# Hemisphere letters
H_ALL = ['N', 'S', 'E', 'W']
H_LAT = ['N', 'S']
H_LON = ['E', 'W']
H_MINUS = ['S', 'W']

def getTmpName():
    """ Creates temporary name based on current time """
    cTime = datetime.datetime.now()               # Current time
    tmpName = str(cTime).replace('-', '')         # Remove hyphens
    tmpName = tmpName.replace(':','')             # Remove colons
    tmpName = tmpName.replace(' ', '_')
    tmpName = tmpName.replace('.', '_')
    tmpName = tmpName[2:] # Trim two first digits of year
    return tmpName
    
def getPolarCoordString(bearing, magVar, distance, distanceUnit, offset, offsetUnit, offsetSide):
    brng = str(bearing)
    MagOrtrue = ''
    if magVar == '' or magVar == 0:
        MagOrTrue = ' TRUE '
    else:
        MagOrTrue = ' MAG '
    dist = str(distance)
    
    if distanceUnit == UOM_M:
        dUnit = ' m'
    elif distanceUnit == UOM_KM:
        dUnit = ' KM'
    elif distanceUnit == UOM_FEET:
        dUnit = ' FEET'
    elif distanceUnit == UOM_SM:
        dUnit = ' SM'
    elif distanceUnit == UOM_NM:
        dUnit = ' NM'
    polarCoordString = brng + MagOrTrue + dist + dUnit + ' ' + offset + ' ' + offsetUnit + ' ' + offsetSide
    return polarCoordString
    
def checkRange(value, valueType):
    """ Check if given value of given value_type is within the range for this value type.
    :param value: float, value to check range
    :valueType : constant of value type (e.g V_AZM)
    """
    if valueType == V_AZM:
        if value < 0 or value > 360:
            result = NOT_VALID
        else:
            result = value
    elif valueType == V_MAG_VAR:
        if value < -360 or value > 360:
            result = NOT_VALID
        else:
            result = value
    elif valueType == V_LAT:
        if (value < -90) or (value > 90):
            result = NOT_VALID
        else:
            result = value
    elif valueType == V_LON:
        if value < -180 or value > 180:
            result = NOT_VALID
        else:
            result = value
    return result
    
# Distance functions
# Patern for distance regular expression
REGEX_DIST = re.compile(r'^\d+(\.\d+)?$') # examples of valid: 0, 0.000, 0.32, 123.455;, examples of invalid: -1.22, s555, 234s5


def checkDistance(d):
    """ Distance validation.
    :param d: string, distance to validate
    :return is_valid: constant VALID if distance is valid, constant NOT_VALID if distance is not valid (e.g distance is less than 0)
    """
    if REGEX_DIST.match(d):
        result = VALID
    else:
        result = NOT_VALID
    return result
    
def km2m(km):
    """ Converts kilometers to meters
    :param km: float, value in kilometers
    :return: value in meters
    """
    return km * 1000

def NM2m(NM):
    """ Converts nautical miles to meters
    :param NM: float, value in nautical miles
    :return: value in meters
    """
    return NM * 1852    

def feet2m(feet):
    """ Converts feet to meters
    :param feet: float, value in feet
    :return: value in meters
    """
    return feet * 0.3048
    
def SM2m(sm):
    """ Converts statue miles to meters
    :param sm: float, value in statue miles
    :return: value in meters
    """
    return sm * 1609.344
    
def toMeters(d, dUnit):
    """ Converts distance given in feet, nautical miles, statue miles etc. to distance in meters
    :param d: float, diatance
    :param dUnit: constant unit of measure, unit of measure
    :return dm: float, distance in meters
    """
    if dUnit == UOM_M:
        dm = d
    elif dUnit == UOM_KM:
        dm = d * 1000
    elif dUnit == UOM_FEET:
        dm = feet2m(d)
    elif dUnit == UOM_SM:
        dm = SM2m(d)
    elif dUnit == UOM_NM:
        dm = NM2m(d)
    return dm



# Azimuth, bearing, magnetic variation functions

def checkAzimuthDD(azm):
    """ Azimuth in DD (decimal degrees) format validation.
    :param azm: string, azimuth to validate
    :return result: float, azimuth value in decima degress format if azm valu is valid azimuth,
                    constant NOT_VALID if azm value is not valid azimuth
    """
    try:
        a = float(azm)
        result = checkRange(a, V_AZM)
    except:
        result = NOT_VALID
    return result

def checkAzimuth(azm):
    """ Azimuth in various format validation.
    :param azm: string, azimuth to validate
    :return result: float, azimuth value in decimal degress format if azm value is valid azimuth,
                    constant NOT_VALID if azm value is not valid azimuth
    """
    azm = str(azm)
    azm.strip()
    azm = azm.replace(",", ".")
    result = checkAzimuthDD(azm)
    return result
    
    
def checkMagVarSignedDD(mag_var):
    """ Magnetic variation  in format signed DD validation.
    Decimal degress format of magnetoc variation, less than 0 - western magnetic variation
    :param mag_var: string, magnetic variation
    :return result: float, value of magbetic variation if  mag_var is correct value or constant NOT_VALID if mag_var
                    is not correct value
    """
    try:
        mv = float(mag_var)
        result = checkRange(mv, V_MAG_VAR)
    except:
        result = NOT_VALID
    return result

def checkMagVarLetterDD(mag_var):
    """ Magnetic variation  in format LDD, DDL validation
    where L depict W (western) or E (eastern) magnetic varaition
    :param mag_var: string, magnetic variation
    :return result: float, value of magbetic variation if  mag_var is correct value or constant NOT_VALID if mag_var
                    is not correct value
    """
    # W15.22, W 15,22, 15.22E, 15.22 EOFError
    if mag_var[0] in ['E', 'W']:
        mag_var_m = mag_var[1:]
        mag_var_m.strip()
        try:
            mv = float(mag_var_m)
            if mag_var_m[0] not in ['0', '1','2', '3', '4', '5', '6', '7', '8', '9']:
                result = NOT_VALID
            elif mv < 0 or mv > 360:  # bnie moze byc np. E-3.5
                result = NOT_VALID
            else:
                if mag_var[0] == 'W':
                    result = -mv
                else:
                    result = mv
        except:
            result = NOT_VALID
        pass
    elif mag_var[len(mag_var) - 1] in ['E', 'W']:
        mag_var_m = mag_var[:len(mag_var) - 2]
        mag_var_m.strip()
        try:
            mv = float(mag_var_m)
            if mag_var_m[0] not in ['0', '1','2', '3', '4', '5', '6', '7', '8', '9']:
                result = NOT_VALID
            elif mv < 0 or mv > 360:  # bnie moze byc np. E-3.5
                result = NOT_VALID
            else:
                if mag_var[0] == 'W':
                    result = -mv
                else:
                    result = mv
        except:
            result = NOT_VALID
    else:
        result = NOT_VALID
    return result

def checkMagVar(mag_var):
    """ Magnetic variation  variosu format validation.
    Decimal degress format of magnetoc variation, less than 0 - western magnetic variation
    :param mag_var: string, magnetic variation
    :return result: float, value of magbetic variation if  mag_var is correct value or constant NOT_VALID if mag_var
                    is not correct value
    """
    mag_var = str(mag_var)
    mag_var.strip()
    mag_var = mag_var.upper()
    mag_var = mag_var.replace(",", ".")
    result = checkMagVarSignedDD(mag_var)
    if result == NOT_VALID:
        result = checkMagVarLetterDD(mag_var)
    return result    
    
# Latitude, longitude functions

def checkSignedDD(dms, cType):
    """ Checks if input parameter is float number, dosen't check latitude, longitude limiest (-90 +90, -180 +180)
    :param dms: string, 
    :param cType: type of coordinate (latitude or longitude)
    return dd: float - decimal degress if is float value, NOT_VALID constant if is not valid float value
    """
    try:
        dd = float(dms)
        dd = checkRange(dd, cType)    
    except:
        dd = NOT_VALID
    return dd
    
def checkHLetterDelimitedDMS_DM(c, cType):
    """ Checks if input parameter is DMS (degrees, minutes, seconds) or DM (degrees, minutes) with hemisphere letter prefix or suffix, 
    :param c: string, coordinate to check
    return dd: float - decimal degress if is valid dms, NOT_VALID constant if is not valid float value
    """
    dd = ''
    pHem = c[0]          
    sHem = c[len(c) - 1] 
    if pHem not in H_ALL and sHem not in H_ALL:
        dd = NOT_VALID
    elif pHem in H_ALL and sHem in H_ALL:
        dd = NOT_VALID
    elif pHem in H_ALL and sHem not in H_ALL:
        dms_n = c[1:] + pHem
    else:
        dms_n = c
    
    if dd != NOT_VALID:
        h = dms_n[len(dms_n) - 1]
        dms_n = dms_n[0:len(dms_n) - 1] # Trim hemisphere letter

        for sep in S_ALL:            # Replace any separator to space separator
             dms_n = dms_n.replace(sep, ' ')
             
        dms_n = dms_n.strip() # Trim trailing spaces - when sconds sign ('') replaced by space at the end of the string
        dms_d = re.sub(r'\s+', ' ', dms_n)  # Replace multiple spaces to single space
        dms_t = dms_d.split(' ')      # Splits dms by spaces and return as tuple
            
        if len(dms_t) == 3:   # 3 elments in tuple - assumes it is DMS format (DD MM SS.sss)
            try:
                d = float(dms_t[0])
                m = float(dms_t[1])
                s = float(dms_t[2])
                if d < 0 or m < 0 or m >= 60 or s < 0 or s >= 60:
                    dd = NOT_VALID
                elif cType == V_LAT and h not in H_LAT:
                    dd = NOT_VALID
                elif cType == V_LON and h not in H_LON:
                    dd = NOT_VALID
                else:
                    dd = d + m/60 + s/3600
                    dd = checkRange(dd, cType) 
                    if (h in H_MINUS) and (dd != NOT_VALID):
                        dd = -dd
            except:
                dd = NOT_VALID
        elif len(dms_t) == 2:   # 2 elments in tuple - assumes it is DM format (DD MM.mmmm)
            try:
                d = float(dms_t[0])
                m = float(dms_t[1])
                if (d < 0) or (m < 0) or (m >= 60):
                    dd = NOT_VALID
                elif cType == V_LAT and h not in H_LAT:
                    dd = NOT_VALID
                elif cType == V_LON and h not in H_LON:
                    dd = NOT_VALID
                else:
                    dd = d + m/60
                    dd = checkRange(dd, cType) 
                    if (h in H_MINUS) and (dd != NOT_VALID):
                        dd = -dd
            except:
                dd = NOT_VALID
        else:
            dd = NOT_VALID
    return dd
    
def parseDMS2DD(dms, cType):
    """ Checks if input parameter is float number, dosen't check latitude, longitude limiest (-90 +90, -180 +180)
    :param dms: string
    :param cType: type of coordinate (latitude or longitude)
    return dd: decimal degrees or NOT_VALID constant if input is not valid coordinate, 
    """
    dms = str(dms)  # Ensure that dms in string variable to perform some string built-in functions
    dms = dms.replace(',','.') # Replace comma decimal separator to period decimal separator
    dms = dms.lstrip(' ')  # Remove leading blanks
    dms = dms.rstrip(' ')  # Remove trailing blanks
    dms = dms.upper()    # Ensure that all charatcers are in upper case, e.g N, E instead of n, e
    
    if dms == '': # Empty string
        dd = NOT_VALID
    else:
        dd = checkSignedDD(dms, cType)
        if dd == NOT_VALID:
            dd = checkHLetterDelimitedDMS_DM(dms, cType)
    return dd
    
def DD2HLetterDelimitedDMS(dd, cType, prec):
    """ Converts coordinate in DD (decimal degrees format) to DMS format space delimited with hemisphere prefix
    :param dd: float, latitude or longitude in decimal degrees format
    :param cType: coordinate type constant, 'LAT' for latitude, 'LON' for longitude
    :param :int, outpur precision (decimal palces of seconds)
    :return dms: string, latitude or longitude in DMS space delimited format with hemisphere letter  - prefix
    """
    
    sec_prec = 3 + prec 
    s_dd = str(dd)
    if s_dd[0] == '-':        # If dd is < 0 remove sign '-'
        s_dd = s_dd[1:]
    
    d = int(math.floor(float(s_dd)))
    m = int(math.floor((float(s_dd) - d) * 60))
    s = (((float(s_dd) - d) * 60) - m) * 60
    
    if cType == V_LAT:
        if d < 10:
            d = '0' + str(d)
        if m < 10:
            m = '0' + str(m)
        if s < 10:
            s = '0' + format(s, '.8f') 
 
        if dd >= 0:
            dms = 'N' + str(d) + ' ' + str(m) + ' ' + str(s)[0:sec_prec]
        else:
            dms = 'S' + str(d) + ' ' + str(m) + ' ' + str(s)[0:sec_prec]
    
    if cType == V_LON:
        if d < 10:
            d = '00' + str(d)
        elif d < 100:
            d = '0' + str(d)
            
        if m < 10:
            m = '0' + str(m)
        if s < 10:
            s = '0' + format(s, '.8f') 
            
        if dd >= 0:
            dms = 'E' + str(d) + ' ' + str(m) + ' ' + str(s)[0:sec_prec]
        else:
            dms = 'W' + str(d) + ' ' + str(m) + ' ' + str(s)[0:sec_prec]
    return dms
    
def vincenty_direct_solution(begin_lat, begin_lon, begin_azimuth, distance, a, b, f): 
    """ Computes the latitude and longitude of the second point based on latitude, longitude,
    of the first point and distance and azimuth from first point to second point.
    Uses the algorithm by Thaddeus Vincenty for direct geodetic problem.
    For more information refer to: http://www.ngs.noaa.gov/PUBS_LIB/inverse.pdf
    
    :param begin_lat: float, latitude of the first point; decimal degrees
    :param begin_lon: float, longitude of the first point; decimal degrees
    :param begin_azimuth: float, azimuth from first point to second point; decimal degrees
    :param distance: float, distance from first point to second point; meters
    :param a: float, semi-major axis of ellispoid; meters
    :param b: float, semi-minor axis of ellipsoid; meters
    :param f: float, flatttening of ellipsoid
    :return lat2_dd, lon2_dd: float, float latitude and longitude of the secon point, decimal degrees
    """
    # Convert latitude, longitude, azimuth of the begining point to radians
    lat1 = math.radians(begin_lat)
    lon1 = math.radians(begin_lon)
    alfa1 = math.radians(begin_azimuth)

    sinAlfa1 = math.sin(alfa1)
    cosAlfa1 = math.cos(alfa1)
    
    # U1 - reduced latitude
    tanU1 = (1 - f) * math.tan(lat1)
    cosU1 = 1 / math.sqrt(1 + tanU1 * tanU1)
    sinU1 = tanU1 * cosU1
    
    # sigma1 - angular distance on the sphere from the equator to begining point
    sigma1 = math.atan2(tanU1, math.cos(alfa1))
    
    # sinAlfa - azimuth of the geodesic at the equator
    sinAlfa = cosU1 * sinAlfa1
    cosSqAlfa = 1 - sinAlfa * sinAlfa
    uSq = cosSqAlfa * (a * a - b * b) / (b * b)
    A = 1 + uSq/16384 * (4096 + uSq * (-768 + uSq * (320 - 175 * uSq)))
    B = uSq/1024 * (256 + uSq * (-128 + uSq * (74 - 47 * uSq)))
    
    sigma = distance / (b * A)
    sigmap = 1
    
    while (math.fabs(sigma - sigmap) > 1e-12):
        cos2sigmaM = math.cos(2 * sigma1 + sigma)
        sinSigma = math.sin(sigma)
        cosSigma = math.cos(sigma)
        dSigma = B*sinSigma*(cos2sigmaM+B/4*(cosSigma*(-1+2*cos2sigmaM*cos2sigmaM)-B/6*cos2sigmaM*(-3+4*sinSigma*sinSigma)*(-3+4*cos2sigmaM*cos2sigmaM)))        
        sigmap = sigma
        sigma = distance / (b * A) + dSigma
    
    var_aux = sinU1 * sinSigma - cosU1 * cosSigma * cosAlfa1 # Auxiliary variable
    
    # Latitude of the end point in radians
    lat2 = math.atan2(sinU1 * cosSigma + cosU1 * sinSigma*cosAlfa1, (1 - f)*math.sqrt(sinAlfa * sinAlfa + var_aux*var_aux))
    
    lamb = math.atan2 (sinSigma * sinAlfa1, cosU1 * cosSigma - sinU1 * sinSigma * cosAlfa1)
    C = f / 16 * cosSqAlfa * (4 + f * (4 - 3 * cosSqAlfa))
    L = lamb - (1 - C) * f * sinAlfa *(sigma + C * sinSigma * (cos2sigmaM + C * cosSigma * (-1 + 2 * cos2sigmaM * cos2sigmaM)))
    # Longitude of the second point in radians
    lon2 = (lon1 + L +3*math.pi)%(2*math.pi) - math.pi
    
    # Convert to decimal degrees
    lat2_dd = math.degrees(lat2)  
    lon2_dd = math.degrees(lon2)
    
    return lat2_dd, lon2_dd
    
def dist_azm_orth_offset2latlon(ref_lat, ref_lon, ref_azm, distance_m, offset_m, offset_side):
    """ Calculates latitude and longitude of the second point base don latitude, longitude of the firts point, azimuth, distance and orthogonal offset
    Example: distance 1500 m, azimuth 45 degress and offset 500 meter left
    :param ref_lat: float, reference point latitude
    :param ref_lon: float, reference poitn longitude
    :param ref_azm: float, azimuth from reference point to intermediate point
    :param distance_m: float, distance in meters
    :param offset_m: float, offset in meters
    :param offset_side: indicate offset side, 'LEFT' for left, 'RIGHT' for right
    :return lat2_dd, lon2_dd: float, second point latitude, longitude
    """
    # Calculate azimuth from intermediate point to second point
    if offset_side == 'LEFT':
        offset_azm = ref_azm - 90
    elif offset_side == 'RIGHT':
        offset_azm = ref_azm + 90
    # Normalize azm to [0,360] degrees
    if offset_azm < 0:
        offset_azm += 360
    elif offset_azm > 360:
        offset_azm -= 360
        
    # Calculate intermediate point latitude, longitude
    inter_lat_dd, inter_lon_dd = vincenty_direct_solution(ref_lat, ref_lon, ref_azm, distance_m, WGS84_A, WGS84_B, WGS84_F)
            
    # Calculate second point latitude, longitude, as reference point use intermediate point
    lat2_dd, lon2_dd = vincenty_direct_solution(inter_lat_dd, inter_lon_dd, offset_azm, offset_m, WGS84_A, WGS84_B, WGS84_F)
    
    return lat2_dd, lon2_dd


def dist_azm_orth_offset2latlon2(ref_lat, ref_lon, x_azm, y_azm, x_m, y_m):
    """ Calculates latitude and longitude of the second point base don latitude, longitude of the firts point, azimuth, distance and orthogonal offset
    Example: distance 1500 m, azimuth 45 degress and offset 500 meter left
    :param ref_lat: float, reference point latitude
    :param ref_lon: float, reference poitn longitude
    :param x_azm: float, X axis azimuth
    :param y_azm: float, Y axis azimuth
    :param x_m: float, distance in meters
    :param y_m: float, offset in meters
    :return lat2_dd, lon2_dd: float, second point latitude, longitude
    """
    # Calculate intermediate point latitude, longitude
    inter_lat_dd, inter_lon_dd = vincenty_direct_solution(ref_lat, ref_lon, x_azm, x_m, WGS84_A, WGS84_B,
                                                          WGS84_F)

    # Calculate second point latitude, longitude, as reference point use intermediate point
    lat2_dd, lon2_dd = vincenty_direct_solution(inter_lat_dd, inter_lon_dd, y_azm, y_m, WGS84_A, WGS84_B,
                                                WGS84_F)

    return lat2_dd, lon2_dd


w = QWidget()


class LatLonByLocCartesian:
    """QGIS Plugin Implementation."""

    def __init__(self, iface):
        """Constructor.

        :param iface: An interface instance that will be passed to this class
            which provides the hook by which you can manipulate the QGIS
            application at run time.
        :type iface: QgisInterface
        """
        self.origLatDD = None
        self.origLonDD = None
        self.origMagVar = None
        self.epName = None
        self.xAxisBearing = None
        self.xAxisAzimuth = None
        self.yAxisAzimuth = None
        self.yAxisOrient = ''
        self.xCoordinate_m = None
        self.yCoordinate_m = None
        self.outLyrSpName = ''
        self.inputFile = ''
        self.outputFile = ''
        self.CSVCoordinateUnit = ''
        
        # Save reference to the QGIS interface
        self.iface = iface
        # initialize plugin directory
        self.plugin_dir = os.path.dirname(__file__)
        # initialize locale
        locale = QSettings().value('locale/userLocale')[0:2]
        locale_path = os.path.join(
            self.plugin_dir,
            'i18n',
            'LatLonByLocCartesian_{}.qm'.format(locale))

        if os.path.exists(locale_path):
            self.translator = QTranslator()
            self.translator.load(locale_path)

            if qVersion() > '4.3.3':
                QCoreApplication.installTranslator(self.translator)

        # Declare instance attributes
        self.actions = []
        self.menu = self.tr(u'&LatLonByLocCartesian')
        # TODO: We are going to let the user set this up in a future iteration
        self.toolbar = self.iface.addToolBar(u'LatLonByLocCartesian')
        self.toolbar.setObjectName(u'LatLonByLocCartesian')

    # noinspection PyMethodMayBeStatic
    def tr(self, message):
        """Get the translation for a string using Qt translation API.

        We implement this ourselves since we do not inherit QObject.

        :param message: String for translation.
        :type message: str, QString

        :returns: Translated version of message.
        :rtype: QString
        """
        # noinspection PyTypeChecker,PyArgumentList,PyCallByClass
        return QCoreApplication.translate('LatLonByLocCartesian', message)


    def add_action(
        self,
        icon_path,
        text,
        callback,
        enabled_flag=True,
        add_to_menu=True,
        add_to_toolbar=True,
        status_tip=None,
        whats_this=None,
        parent=None):
        """Add a toolbar icon to the toolbar.

        :param icon_path: Path to the icon for this action. Can be a resource
            path (e.g. ':/plugins/foo/bar.png') or a normal file system path.
        :type icon_path: str

        :param text: Text that should be shown in menu items for this action.
        :type text: str

        :param callback: Function to be called when the action is triggered.
        :type callback: function

        :param enabled_flag: A flag indicating if the action should be enabled
            by default. Defaults to True.
        :type enabled_flag: bool

        :param add_to_menu: Flag indicating whether the action should also
            be added to the menu. Defaults to True.
        :type add_to_menu: bool

        :param add_to_toolbar: Flag indicating whether the action should also
            be added to the toolbar. Defaults to True.
        :type add_to_toolbar: bool

        :param status_tip: Optional text to show in a popup when mouse pointer
            hovers over the action.
        :type status_tip: str

        :param parent: Parent widget for the new action. Defaults None.
        :type parent: QWidget

        :param whats_this: Optional text to show in the status bar when the
            mouse pointer hovers over the action.

        :returns: The action that was created. Note that the action is also
            added to self.actions list.
        :rtype: QAction
        """

        # Create the dialog (after translation) and keep reference
        self.dlg = LatLonByLocCartesianDialog()

        icon = QIcon(icon_path)
        action = QAction(icon, text, parent)
        action.triggered.connect(callback)
        action.setEnabled(enabled_flag)

        if status_tip is not None:
            action.setStatusTip(status_tip)

        if whats_this is not None:
            action.setWhatsThis(whats_this)

        if add_to_toolbar:
            self.toolbar.addAction(action)

        if add_to_menu:
            self.iface.addPluginToMenu(
                self.menu,
                action)

        self.actions.append(action)

        return action

    def initGui(self):
        """Create the menu entries and toolbar icons inside the QGIS GUI."""

        icon_path = ':/plugins/LatLonByLocCartesian/icon.png'
        self.add_action(
            icon_path,
            text=self.tr(u'LatLonByLocCartesian'),
            callback=self.run,
            parent=self.iface.mainWindow())

        self.dlg.pbAddPoint.clicked.connect(self.calcSpPoint)
        self.dlg.pbInCsv.clicked.connect(self.selectInputFile)
        self.dlg.pbOutCsv.clicked.connect(self.selectOutputFile)
        self.dlg.pbCalcCSV.clicked.connect(self.calcCSVPoints)

    def unload(self):
        """Removes the plugin menu item and icon from QGIS GUI."""
        for action in self.actions:
            self.iface.removePluginMenu(
                self.tr(u'&LatLonByLocCartesian'),
                action)
            self.iface.removeToolBarIcon(action)
        # remove the toolbar
        del self.toolbar
    
    def getYAxisOrientation(self):
        """ Gets y axis orientation """
        if self.dlg.cbyAxisOrient.currentText() == 'left':
            side = 'LEFT'
        elif self.dlg.yAxisOrient.currentText() == 'right':
            side = 'RIGHT'
        return side

    def createTmpLayer(self, lyrName):
        """ Create temporary 'memory' layer to store results of calculations
        :param lyrName: string, layer name
        """
        tmpLyr = QgsVectorLayer('Point?crs=epsg:4326', lyrName, 'memory')
        prov = tmpLyr.dataProvider()
        tmpLyr.startEditing()
        prov.addAttributes([QgsField("PNAME", QVariant.String),
                            QgsField("LAT_DMS", QVariant.String),
                            QgsField("LON_DMS", QVariant.String),
                            QgsField("POLAR_COOR", QVariant.String)])
        tmpLyr.commitChanges()
        QgsMapLayerRegistry.instance().addMapLayers([tmpLyr])

    """ Single point functions """
    
    def getSpXUnit(self):
        """ Gets X coordinate unit """
        if self.dlg.cbSpXUnit.currentText() == 'm':
            unit = UOM_M
        elif self.dlg.cbSpXUnit.currentText() == 'KM':
            unit = UOM_KM
        elif self.dlg.cbSpXUnit.currentText() == 'NM':
            unit = UOM_NM
        elif self.dlg.cbSpXUnit.currentText() == 'feet':
            unit = UOM_FEET
        elif self.dlg.cbSpXUnit.currentText() == 'SM':
            unit = UOM_SM
        return unit

    def getSpYUnit(self):
        """ Gets Y coordinate unit """
        if self.dlg.cbSpYUnit.currentText() == 'm':
            unit = UOM_M
        elif self.dlg.cbSpYUnit.currentText() == 'KM':
            unit = UOM_KM
        elif self.dlg.cbSpYUnit.currentText() == 'NM':
            unit = UOM_NM
        elif self.dlg.cbSpYUnit.currentText() == 'feet':
            unit = UOM_FEET
        elif self.dlg.cbSpYUnit.currentText() == 'SM':
            unit = UOM_SM
        return unit    
        
    def checkSpLocCartesianInput(self):
        checkResult = True
        errMsg = ''
        
        origLatDMS = self.dlg.leOrigLat.text()        # Latitude of the reference point
        origLonDMS = self.dlg.leOrigLon.text()        # Longitude of the reference point
        magVar = self.dlg.leOrigMagVar.text()    # Magnetic variation of the reference point
        self.epName = self.dlg.leEpName.text() 
        self.xAxisBearing = self.dlg.lexAxisBearing.text()
        self.yAxisOrient = self.getYAxisOrientation()
        x = float(self.dlg.leEpX.text())
        y = float(self.dlg.leEpY.text())

        if origLonDMS == '':
            errMsg += 'Enter origin latitude\n'
        else:
            self.origLatDD = parseDMS2DD(origLatDMS, V_LAT)
            if self.origLatDD == NOT_VALID:
                errMsg += 'Latitude not valid\n'
                checkResult = False

        if origLonDMS == '':
            errMsg += 'Enter origin longitude\n'
        else:
            self.origLonDD = parseDMS2DD(origLonDMS, V_LON)
            if self.origLatDD == NOT_VALID:
                errMsg += 'Longitude not valid\n'
                checkResult = False

        if magVar == '':  # Magnetic Variation not entered - assume magnetic variation as 0.0
            self.origMagVar = 0.0
        else:
            self.origMagVar = checkMagVar(magVar)
            if self.origMagVar == NOT_VALID:
                errMsg += 'Magnetic variation wrong format!\n'
                checkResult = False

        if checkResult:
            self.xAxisAzimuth = float(self.xAxisBearing) + self.origMagVar
            if self.xAxisAzimuth < 0:
                self.xAxisAzimuth += 360
            elif self.xAxisAzimuth > 360:
                self.xAxisAzimuth -= 360

            self.yAxisOrient = self.getYAxisOrientation()
            if self.yAxisOrient == 'RIGHT':
                self.yAxisAzimuth = self.xAxisAzimuth + 90
            elif self.yAxisOrient == 'LEFT':
                self.yAxisAzimuth = self.xAxisAzimuth - 90

            if self.yAxisAzimuth < 0:
                self.yAxisAzimuth += 360
            elif self.yAxisAzimuth > 360:
                self.yAxisAzimuth -= 360

            xUnit = self.getSpXUnit()
            yUnit = self.getSpYUnit()

            self.xCoordinate_m = toMeters(x, xUnit)
            self.yCoordinate_m = toMeters(y, yUnit)
        else:
            QMessageBox.critical(w, "Message", errMsg)

        return checkResult
    
    def calcSpPoint(self):
        if self.checkSpLocCartesianInput():
            xAxisAzimuthReverse = self.xAxisAzimuth - 180
            if xAxisAzimuthReverse < 0:
                xAxisAzimuthReverse += 360
            elif xAxisAzimuthReverse > 360:
                xAxisAzimuthReverse -= 360

            yAxisAzimuthReverse = self.yAxisAzimuth - 180
            if yAxisAzimuthReverse < 0:
                yAxisAzimuthReverse += 360
            elif yAxisAzimuthReverse > 360:
                yAxisAzimuthReverse -= 360


            if self.xCoordinate_m >= 0: # Normal x axis azimuth
                if self.yCoordinate_m >= 0: # Normal y axis azimuth
                    epLatDD, epLonDD = dist_azm_orth_offset2latlon2(self.origLatDD, self.origLonDD,
                                                                    self.xAxisAzimuth, self.yAxisAzimuth,
                                                                    self.xCoordinate_m, self.yCoordinate_m)
                elif self.yCoordinate_m < 0:
                    epLatDD, epLonDD = dist_azm_orth_offset2latlon2(self.origLatDD, self.origLonDD,
                                                                    self.xAxisAzimuth, yAxisAzimuthReverse,
                                                                    self.xCoordinate_m, math.fabs(self.yCoordinate_m))
            elif self.xCoordinate_m < 0:
                if self.yCoordinate_m >= 0: # Normal y axis azimuth
                    epLatDD, epLonDD = dist_azm_orth_offset2latlon2(self.origLatDD, self.origLonDD,
                                                                    xAxisAzimuthReverse, self.yAxisAzimuth,
                                                                    math.fabs(self.xCoordinate_m), self.yCoordinate_m)
                elif self.yCoordinate_m < 0:
                    epLatDD, epLonDD = dist_azm_orth_offset2latlon2(self.origLatDD, self.origLonDD,
                                                                    xAxisAzimuthReverse, yAxisAzimuthReverse,
                                                                    math.fabs(self.xCoordinate_m), math.fabs(self.yCoordinate_m))
            epLatDMS = DD2HLetterDelimitedDMS(epLatDD, V_LAT, 3)
            epLonDMS = DD2HLetterDelimitedDMS(epLonDD, V_LON, 3)
            polarCoord = 'TEST'

            layers = self.iface.legendInterface().layers()
            layerList = []  # List of layers in current (opened) QGIS project
            for layer in layers:
                layerList.append(layer.name())

            if self.outLyrSpName not in layerList:
                self.outLyrSpName = getTmpName()
                self.createTmpLayer(self.outLyrSpName)
                outLyr = QgsVectorLayer('Point?crs=epsg:4326', self.outLyrSpName, 'memory')
                outLyr = self.iface.activeLayer()
                outLyr.startEditing()
                outProv = outLyr.dataProvider()
                feat = QgsFeature()
                # Add reference point to layer
                feat.setGeometry(QgsGeometry.fromPoint(QgsPoint(self.origLonDD, self.origLatDD)))
                feat.setAttributes(['REF_POINT', self.dlg.leOrigLat.text(), self.dlg.leOrigLon.text()])
                outProv.addFeatures([feat])
                outLyr.commitChanges()
                # Add calculated point to layer
                feat.setGeometry(QgsGeometry.fromPoint(QgsPoint(epLonDD, epLatDD)))
                feat.setAttributes([self.epName, epLatDMS, epLonDMS, polarCoord])
                outProv.addFeatures([feat])
                outLyr.commitChanges()
                outLyr.updateExtents()
                self.iface.mapCanvas().setExtent(outLyr.extent())
                self.iface.mapCanvas().refresh()
            elif self.outLyrSpName in layerList:
                outLyr = QgsVectorLayer('Point?crs=epsg:4326', self.outLyrSpName, 'memory')
                outLyr = self.iface.activeLayer()
                outLyr.startEditing()
                outProv = outLyr.dataProvider()
                feat = QgsFeature()
                feat.setGeometry(QgsGeometry.fromPoint(QgsPoint(epLonDD, epLatDD)))
                feat.setAttributes([self.epName, epLatDMS, epLonDMS, polarCoord])
                outProv.addFeatures([feat])
                outLyr.commitChanges()
                outLyr.updateExtents()
                self.iface.mapCanvas().setExtent(outLyr.extent())
                self.iface.mapCanvas().refresh()
        return

    def selectInputFile(self):
        """ Select input csv file with data: ID of the point, azimuth, distance from refernce point to the current point """
        self.inputFile = QFileDialog.getOpenFileName(self.dlg, "Select input file ", "", '*.csv')
        self.dlg.leInCSV.setText(self.inputFile)

    def selectOutputFile(self):
        """ Select output csv file """
        self.outputFile = QFileDialog.getSaveFileName(self.dlg, "Select output file ", "", '*.csv')
        self.dlg.leOutCSV.setText(self.outputFile)

    def getCSVCoordinateUnit(self):
        """ Gets X, Y coordinates unit for CSV file data """
        if self.dlg.cbCSVCoorUnit.currentText() == 'm':
            unit = UOM_M
        elif self.dlg.cbCSVCoorUnit.currentText() == 'KM':
            unit = UOM_KM
        elif self.dlg.cbCSVCoorUnit.currentText() == 'NM':
            unit = UOM_NM
        elif self.dlg.cbCSVCoorUnit.currentText() == 'feet':
            unit = UOM_FEET
        elif self.dlg.cbCSVCoorUnit.currentText() == 'SM':
            unit = UOM_SM
        return unit

    def checkCSVPoints(self):
        checkResult = True
        errMsg = ''

        origLatDMS = self.dlg.leOrigLat.text()  # Latitude of the reference point
        origLonDMS = self.dlg.leOrigLon.text()  # Longitude of the reference point
        magVar = self.dlg.leOrigMagVar.text()  # Magnetic variation of the reference point
        self.xAxisBearing = self.dlg.lexAxisBearing.text()
        self.yAxisOrient = self.getYAxisOrientation()

        self.inputFile = self.dlg.leInCSV.text()
        self.outputFile = self.dlg.leOutCSV.text()


        if origLonDMS == '':
            errMsg += 'Enter origin latitude\n'
        else:
            self.origLatDD = parseDMS2DD(origLatDMS, V_LAT)
            if self.origLatDD == NOT_VALID:
                errMsg += 'Latitude not valid\n'
                checkResult = False

        if origLonDMS == '':
            errMsg += 'Enter origin longitude\n'
        else:
            self.origLonDD = parseDMS2DD(origLonDMS, V_LON)
            if self.origLatDD == NOT_VALID:
                errMsg += 'Longitude not valid\n'
                checkResult = False

        if magVar == '':  # Magnetic Variation not entered - assume magnetic variation as 0.0
            self.origMagVar = 0.0
        else:
            self.origMagVar = checkMagVar(magVar)
            if self.origMagVar == NOT_VALID:
                errMsg += 'Magnetic variation wrong format!\n'
                checkResult = False

        if self.inputFile == '':
            errMsg += 'Choose input file\n'
            checkResult = False
        # TO DO: check input file header

        if self.outputFile == '':
            errMsg += 'Choose output file\n'
            checkResult = False

        self.CSVCoordinateUnit = self.getCSVCoordinateUnit()

        if checkResult:
            self.xAxisAzimuth = float(self.xAxisBearing) + self.origMagVar
            if self.xAxisAzimuth < 0:
                self.xAxisAzimuth += 360
            elif self.xAxisAzimuth > 360:
                self.xAxisAzimuth -= 360

            self.yAxisOrient = self.getYAxisOrientation()
            if self.yAxisOrient == 'RIGHT':
                self.yAxisAzimuth = self.xAxisAzimuth + 90
            elif self.yAxisOrient == 'LEFT':
                self.yAxisAzimuth = self.xAxisAzimuth - 90

            if self.yAxisAzimuth < 0:
                self.yAxisAzimuth += 360
            elif self.yAxisAzimuth > 360:
                self.yAxisAzimuth -= 360
        else:
            QMessageBox.critical(w, "Message", errMsg)
        return checkResult

    def calcCSVPoints(self):
        if self.checkCSVPoints():
            errNotes = ''
            outLyrCSVName = 'CSV' + getTmpName()
            self.createTmpLayer(outLyrCSVName)
            outLyrCSV = QgsVectorLayer('Point?crs=epsg:4326', outLyrCSVName, 'memory')
            outLyrCSV = self.iface.activeLayer()
            outLyrCSV.startEditing()
            outProvCSV = outLyrCSV.dataProvider()

            feat = QgsFeature()

            outCSVFieldNames = ['P_NAME',
                                'X',
                                'Y',
                                'LAT_DMS',
                                'LON_DMS',
                                'CART_COOR_STRING',
                                'ERR_NOTES']

            xAxisAzimuthReverse = self.xAxisAzimuth - 180
            if xAxisAzimuthReverse < 0:
                xAxisAzimuthReverse += 360
            elif xAxisAzimuthReverse > 360:
                xAxisAzimuthReverse -= 360

            yAxisAzimuthReverse = self.yAxisAzimuth - 180
            if yAxisAzimuthReverse < 0:
                yAxisAzimuthReverse += 360
            elif yAxisAzimuthReverse > 360:
                yAxisAzimuthReverse -= 360

            with open(self.inputFile, 'r') as inCSV:
                with open(self.outputFile, 'w') as outCSV:
                    reader = csv.DictReader(inCSV, delimiter=';')
                    writer = csv.DictWriter(outCSV, fieldnames=outCSVFieldNames, delimiter=';')
                    for row in reader:
                        try:
                            # TO DO: validation input line
                            x_m = toMeters(float(row['X']), self.CSVCoordinateUnit)
                            y_m = toMeters(float(row['Y']), self.CSVCoordinateUnit)

                            if x_m >= 0:  # Normal x axis azimuth
                                if y_m >= 0:  # Normal y axis azimuth
                                    epLatDD, epLonDD = dist_azm_orth_offset2latlon2(self.origLatDD, self.origLonDD,
                                                                                    self.xAxisAzimuth,
                                                                                    self.yAxisAzimuth,
                                                                                    x_m, y_m)
                                elif y_m < 0:
                                    epLatDD, epLonDD = dist_azm_orth_offset2latlon2(self.origLatDD, self.origLonDD,
                                                                                    self.xAxisAzimuth,
                                                                                    yAxisAzimuthReverse,
                                                                                    x_m, math.fabs(y_m))
                            elif x_m < 0:
                                if y_m >= 0:  # Normal y axis azimuth
                                    epLatDD, epLonDD = dist_azm_orth_offset2latlon2(self.origLatDD, self.origLonDD,
                                                                                    xAxisAzimuthReverse,
                                                                                    self.yAxisAzimuth,
                                                                                    math.fabs(x_m), y_m)
                                elif y_m < 0:
                                    epLatDD, epLonDD = dist_azm_orth_offset2latlon2(self.origLatDD, self.origLonDD,
                                                                                    xAxisAzimuthReverse,
                                                                                    yAxisAzimuthReverse,
                                                                                    math.fabs(self.x_m),
                                                                                    math.fabs(self.y_m))
                            epLatDMS = DD2HLetterDelimitedDMS(epLatDD, V_LAT, 3)
                            epLonDMS = DD2HLetterDelimitedDMS(epLonDD, V_LON, 3)

                            cartCoor = 'TEST'

                            writer.writerow({'P_NAME': row['P_NAME'],
                                             'X': row['X'],
                                             'Y': row['Y'],
                                             'LAT_DMS': epLatDMS,
                                             'LON_DMS': epLonDMS,
                                             'CART_COOR_STRING': 'TEST',
                                             'ERR_NOTES': ''})
                            endPoint = QgsPoint(epLonDD, epLatDD)
                            feat.setGeometry(QgsGeometry.fromPoint(endPoint))
                            feat.setAttributes([row['P_NAME'], epLatDMS, epLonDMS, cartCoor])
                            outProvCSV.addFeatures([feat])
                            outLyrCSV.commitChanges()
                        except:
                            pass
            outLyrCSV.updateExtents()
            self.iface.mapCanvas().setExtent(outLyrCSV.extent())
            self.iface.mapCanvas().refresh()
        return

    def run(self):
        """Run method that performs all the real work"""
        # show the dialog
        self.dlg.show()
        # Run the dialog event loop
        result = self.dlg.exec_()
        # See if OK was pressed
        if result:
            # Do something useful here - delete the line containing pass and
            # substitute with your code.
            pass

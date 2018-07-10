# -*- coding: utf-8 -*-
"""
/***************************************************************************
 LatLonByLocCartesian
                                 A QGIS plugin
 latitude, longitude calcucation based on local cartesian coordinate system
                             -------------------
        begin                : 2018-07-09
        copyright            : (C) 2018 by Paweł Strzelewicz
        email                : pawel.strzelewicz@wp.pl
        git sha              : $Format:%H$
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
 This script initializes the plugin, making it known to QGIS.
"""


# noinspection PyPep8Naming
def classFactory(iface):  # pylint: disable=invalid-name
    """Load LatLonByLocCartesian class from file LatLonByLocCartesian.

    :param iface: A QGIS interface instance.
    :type iface: QgisInterface
    """
    #
    from .latlon_by_loc_cartesian import LatLonByLocCartesian
    return LatLonByLocCartesian(iface)

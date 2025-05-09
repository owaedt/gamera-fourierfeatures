#!/usr/bin/env python
#
# Copyright (C) 2013 Christian Brandt
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#  
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#

from pathlib import Path

from setuptools import setup

from gamera import gamera_setup

# This constant should be the name of the toolkit
TOOLKIT_NAME = "fourierfeatures"

# ----------------------------------------------------------------------------
# You should not usually have to edit anything below, but it is
# implemented here and not in the Gamera core so that you can edit it
# if you need to do something more complicated (for example, building
# and linking to a third- party library).
# ----------------------------------------------------------------------------

PLUGIN_PATH = 'gamera/toolkits/%s/plugins/' % TOOLKIT_NAME
PACKAGE = 'gamera.toolkits.%s' % TOOLKIT_NAME
PLUGIN_PACKAGE = PACKAGE + ".plugins"
plugins = gamera_setup.get_plugin_filenames(PLUGIN_PATH)
plugin_extensions = gamera_setup.generate_plugins(plugins, PLUGIN_PACKAGE)

# This is a standard setuptools setup initializer. If you need to do
# anything more complex here, refer to the Python setuptools documentation.
if __name__ == "__main__":
    setup(name="gamera-fourierfeatures",
          version="2.0.0",
          ext_modules=plugin_extensions,
          packages=[PACKAGE, PLUGIN_PACKAGE],
          include_dirs=['include/plugins'],
          python_requires='>=3.5',
          scripts=[],
          install_requires=['gamera>=4.1.0']
          )

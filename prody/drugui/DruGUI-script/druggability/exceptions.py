# Druggability: Python Package and VMD GUI for Druggability Index Analysis
# Copyright (C) 2010  Ahmet Bakan <ahb12@pitt.edu>
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010  Ahmet Bakan'
__version__ = '0.5.2'


"""Druggability package related exceptions."""

class GridError(Exception):    

    """Exception for errors related to handling grid files and data."""


    pass
    
    
class ProbeError(Exception):    

    """Exception for errors related to handling probe grid and data."""


    pass
    
    
class DIAError(Exception):

    """Exception for errors related to druggability index calculations."""


    pass

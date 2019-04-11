# -*- coding: utf-8 -*-
"""
Exceptions.
"""

class XRRError(Exception):
    """
    XRR base exception.
    """
    pass
    
    
class XRRSystemError(XRRError):
    """
    Base class for system errors.
    """
    pass
    
    
class XRRLibraryVersionError(XRRSystemError):
    """
    Raised whenever an external library does not have the correct version.
    """
    def __init__(self, library, min_version, current_version):
        super(XRRLibraryVersionError, self).__init__()
        self._lib = library
        self._min_v = min_version
        self._cur_v = current_version
        
    def __str__(self):
        msg = "Version {0} of '{1}' required. You are using version {2}."
        return msg.format(self._min_v, self._lib, self._cur_v)

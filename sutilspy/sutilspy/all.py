#!/usr/bin/env python
# Copyright (C) 2017 Sur Herrera Paredes

# import subprocess
# import os

class Error(Exception):
    """Base class for exceptions in this module."""
    pass

class IntegrityError(Error):
    """Exception raised for failure in the vdb-vvalidate."""
    pass

class MissingFileError(Error):
    """Exception raised for missing files."""
    pass

class ProcessError(Error):
    """Exception raised for failure in the some sra-tools call."""
    pass




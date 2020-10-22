"""AMICI-generated module for model Leber2015"""

import amici

# Ensure we are binary-compatible, see #556
if '0.10.7' != amici.__version__:
    raise RuntimeError('Cannot use model Leber2015, generated with AMICI '
                       'version 0.10.7, together with AMICI version'
                       f' {amici.__version__} which is present in your '
                       'PYTHONPATH. Install the AMICI package matching the '
                       'model version or regenerate the model with the AMICI '
                       'currently in your path.')

from Leber2015.Leber2015 import *

__version__ = '0.1.0'

"""AMICI-generated module for model model2_"""

import amici

# Ensure we are binary-compatible, see #556
if '0.10.7' != amici.__version__:
    raise RuntimeError('Cannot use model model2_, generated with AMICI '
                       'version 0.10.7, together with AMICI version'
                       f' {amici.__version__} which is present in your '
                       'PYTHONPATH. Install the AMICI package matching the '
                       'model version or regenerate the model with the AMICI '
                       'currently in your path.')

from model2_.model2_ import *

__version__ = '0.1.0'

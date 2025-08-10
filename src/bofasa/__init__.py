"""
Bacterial Orthology Finding And Syntenic Analysis (bofasa)

A tool for high-quality orthology inference across multiple bacterial species.
"""

__version__ = "1.2.0"
__author__ = "Rauf Salamzade"
__email__ = "salamzader@gmail.com"

# Import only lightweight modules by default
from .config import *

# Import heavy modules only when needed
# from .core import *
# from .utils import *
# from .analysis import *
# from .processing import *

__all__ = [
    # Configuration
    "get_version",
    "reset_global_variables",
]

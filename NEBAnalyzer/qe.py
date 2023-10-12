import os
import numpy as np
from ase import io

from . import Analyzer


class QEAnalyzer(Analyzer):

    def __init__(self, ddir: str, read_neb: bool = True):

        super().__init__(ddir, read_neb)

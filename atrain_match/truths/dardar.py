# -*- coding: utf-8 -*-
# Copyright (c) 2009-2019 atrain_match developers
#
# This file is part of atrain_match.
#
# atrain_match is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# atrain_match is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with atrain_match.  If not, see <http://www.gnu.org/licenses/>.
# Change log is found in git
from atrain_match.utils.runutils import do_some_logging
from atrain_match.truths.calipso import (find_break_points)
from atrain_match.libs.extract_imager_along_track import imager_track_from_matched
from atrain_match.utils.common import (MatchupError, ProcessingError,
                                       elements_within_range)
import atrain_match.config as config
from atrain_match.matchobject_io import (CloudsatObject,
                                         TruthImagerTrackObject)
import time
import numpy as np
import os
import logging
logger = logging.getLogger(__name__)

def get_dardar(filename):
    pass

def read_dardar_nc(filename):
    pass
#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) 2019 atrain_match developers
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
# Author(s):

#   Nina Hakansson <nina.hakansson@smhi.se>

"""Test Initializer for atrain_match."""

import unittest

from atrain_match.tests import (test_add_single_shot,
                                test_detection_height,
                                test_init,
                                test_plot_kuipers_on_map,
                                test_utils)


def suite():
    """Test global test suite."""
    mysuite = unittest.TestSuite()
    mysuite.addTests(test_add_single_shot.suite())
    mysuite.addTests(test_detection_height.suite())
    mysuite.addTests(test_init.suite())
    mysuite.addTests(test_plot_kuipers_on_map.suite())
    mysuite.addTests(test_utils.suite())
    return mysuite


def load_tests(loader, tests, pattern):
    """Load all tests."""
    return suite()

#!/usr/bin/env python3

import os, sys, platform
import multiprocessing as mp
import pyctest.pyctest as pyctest
import pyctest.helpers as helpers

directory = "./"

# these are required

pyctest.PROJECT_NAME = "AORSA"
pyctest.SOURCE_DIRECTORY = directory
pyctest.BINARY_DIRECTORY = directory

args = helpers.ArgumentParser(pyctest.PROJECT_NAME,
                              pyctest.SOURCE_DIRECTORY,
                              pyctest.BINARY_DIRECTORY).parse_args()

pyctest.BUILD_COMMAND = "make"

test = pyctest.test()
test.SetName("Helicon")
test.SetProperty("WORKING_DIRECTORY","test/DIIID-helicon")
test.SetCommand(["./test.sh"])

test = pyctest.test()
test.SetName("Whistler")
test.SetProperty("WORKING_DIRECTORY","test/DIIID-whistler")
test.SetCommand(["./test.sh"])

pyctest.run()

#!/usr/bin/env python

import os, sys, platform
import multiprocessing as mp
import pyctest.pyctest as pyctest
import pyctest.helpers as helpers

parser = helpers.ArgumentParser("AORSA", source_dir=os.getcwd(), binary_dir=os.getcwd(),
                                vcs_type="git")
parser.add_argument("-n", "--build", type=str, required=True, help="Build name for identification")
args = parser.parse_args()

# CONFIGURE_COMMAND can only run one command so if autogen is required, just execute it here
#cmd = pyctest.command(["./autogen.sh"])
#cmd.SetWorkingDirectory(pyctest.SOURCE_DIRECTORY)
#cmd.SetErrorQuiet(False)
#cmd.Execute()

pyctest.BUILD_NAME = "{}".format(args.build)
#pyctest.CONFIGURE_COMMAND = "./configure"
pyctest.BUILD_COMMAND = "make"

test = pyctest.test()
test.SetName("Helicon")
test.SetProperty("WORKING_DIRECTORY",pyctest.SOURCE_DIRECTORY+"/test/DIIID-helicon")
test.SetCommand(["./test.sh"])

test = pyctest.test()
test.SetName("Whistler")
test.SetProperty("WORKING_DIRECTORY",pyctest.SOURCE_DIRECTORY+"/test/DIIID-whistler")
test.SetCommand(["./test.sh"])

pyctest.run()

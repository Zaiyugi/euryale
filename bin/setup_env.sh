#!/bin/bash

# Post-make script for editing environment variables

# Create an absolute path to the project directory
ABSOLUTE_PATH=$(cd .; pwd)

# Check for required libraries in path
# If not found, add to LD_LIBRARY_PATH
if [[ $LD_LIBRARY_PATH != *"group/dpa/lib"* ]]
then
	export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ABSOLUTE_PATH:/group/dpa/lib
fi


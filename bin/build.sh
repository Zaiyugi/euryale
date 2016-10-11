#!/bin/bash

# Build script for Cauldron repo needed by project

# Create an absolute path to the project directory
ABSOLUTE_PATH=$(cd .; pwd)

# Build everything
echo "*** Building required libraries... ***"
echo "--------------------------------------"
make all
echo $'-------------------------------------\n'

# Build project
echo "*** Building project: $1 ***"
echo "----------------------------"
#echo $1
make $1
echo $'---------------------------\n'

# Add necessary environment variables
echo "*** Setting up environment ***"
. ./bin/setup_env.sh


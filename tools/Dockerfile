# This Dockerfile creates the polymec/polymec-dev:latest image on Docker Hub.
# It was created and uploaded with the following commands (with docker running):
#
# docker build -t polymec/polymec-dev:latest .
# docker login
# docker push polymec/polymec-dev:latest

# Use ubuntu:bionic as a base image
FROM ubuntu:bionic

# Credit/blame where it's due
LABEL maintainer="pbtoast"

# Update the apt repository and install what we need.
RUN apt-get update
RUN apt-get install -y openmpi-bin libopenmpi-dev ssh g++ liblapack-dev cmake ninja-build autotools-dev bash git valgrind libpetsc3.7.7 lcov curl;

# Symlink petsc so we can use it easily with Travis CI.
RUN ln -s /usr/lib/x86_64-linux-gnu/libpetsc_real.so.3.7.7 /usr/lib/libpetsc.so

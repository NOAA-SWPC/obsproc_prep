#!/bin/sh
set -x
set -e    # fail if an error is hit so that errors do not go unnoticed
if [ $# -eq 0 ]; then
  dir_list=*.fd
else
  dir_list=$*
fi
echo $dir_list
#source ./load_libs.rc  # use modules to set library related environment variables
source ./setlibs.rc  # use this if existing library modules don't quite cover all that is needed.

for sdir in $dir_list; do
 dir=${sdir%\/}  # chop trailing slash if necessary
 [ $dir = bufr_tranwindsat.fd ] && continue  # this code isn't ready
 cd $dir
 make clobber
#  uncomment following two lines if using modules (load_libs.rc)
# [ $dir = bufr_trantmi.fd ] && module load HDF4/$HDF4_ver 
# [ $dir = bufr_tranmls.fd -o $dir = bufr_tranomi.fd ] && module load HDF5/$HDF5_ver/serial
 make
#  uncomment following two lines if using modules (load_libs.rc)
# [ $dir = bufr_trantmi.fd ] && module unload HDF4/$HDF4_ver
# [ $dir = bufr_tranmls.fd -o $dir = bufr_tranomi.fd ] && module unload HDF5/$HDF5_ver/serial
 cd ..
done



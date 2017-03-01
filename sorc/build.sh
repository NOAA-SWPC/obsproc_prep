#!/bin/sh
set -x
set -e    # fail if an error is hit so that errors do not go unnoticed

##  determine system/phase   (update for theia later -dcs)
module load prod_util
sys_tp=$(getsystem.pl -tp)
echo $sys_tp

module purge
case $sys_tp in
 IBM-p1|IBM-p2)
   module load ics/12.1;
   module load ibmpe;
   ;;
 Cray-XC40)
   module load PrgEnv-intel;
   module load craype-haswell;
   module load cray-mpich/7.2.0;
   module swap intel/16.3.210;
   module load iobuf/2.0.7;
   lib_build="intel";
   export FC=ftn;
   ;;
 *) echo unexpected system.  Update for $sys_tp;;
esac

source ./load_libs.rc  # use modules to set library related environment variables
#source ./setlibs.rc  # use this if existing library modules don't quite cover all that is needed.

module list

if [ $# -eq 0 ]; then
  dir_list=*.fd
else
  dir_list=$*
fi
echo $dir_list

clobber=${clobber:-clobber_yes}  # user can override the default of running "make clobber"
for sdir in $dir_list; do
 dir=${sdir%\/}  # chop trailing slash if necessary
 cd $dir
 [ $clobber != clobber_no ]  && make clobber
 if [ $sys_tp = Cray-XC40 ]; then
   make FC=$FC
 else
   make
 fi
 ###touch *
 ls -l
 cd ..
done



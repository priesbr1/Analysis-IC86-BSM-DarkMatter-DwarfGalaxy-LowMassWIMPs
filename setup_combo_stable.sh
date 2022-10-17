#!/bin/bash

eval `/cvmfs/icecube.opensciencegrid.org/py3-v4.1.0/setup.sh`
module purge

export PATH=/opt/software/core/lua/lua/bin:/cvmfs/icecube.opensciencegrid.org/py3-v4.1.0/RHEL_7_x86_64/metaprojects/combo/stable/bin:/cvmfs/icecube.opensciencegrid.org/py3-v4.1.0/RHEL_7_x86_64/bin:/usr/lib64/qt-3.3/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/usr/local/hpcc/bin:/usr/lpp/mmfs/bin:/opt/ibutils/bin:/opt/puppetlabs/bin
export LD_LIBRARY_PATH=/cvmfs/icecube.opensciencegrid.org/py3-v4.1.0/RHEL_7_x86_64/metaprojects/combo/stable/lib:/cvmfs/icecube.opensciencegrid.org/py3-v4.1.0/RHEL_7_x86_64/metaprojects/combo/stable/cernroot/lib:/cvmfs/icecube.opensciencegrid.org/py3-v4.1.0/RHEL_7_x86_64/metaprojects/combo/stable/lib/tools:/cvmfs/icecube.opensciencegrid.org/py3-v4.1.0/RHEL_7_x86_64/spack/opt/spack/linux-centos7-x86_64/gcc-9.2.0spack/geant4-10.04-jtcfrttq4krvbvabfjtkq7itsf3dieqg/lib64:/cvmfs/icecube.opensciencegrid.org/py3-v4.1.0/RHEL_7_x86_64/lib:/cvmfs/icecube.opensciencegrid.org/py3-v4.1.0/RHEL_7_x86_64/lib64:/tmp/tmpkvezwnnt/lib:/usr/lib64/nvidia
export DYLD_LIBRARY_PATH=/cvmfs/icecube.opensciencegrid.org/py3-v4.1.0/RHEL_7_x86_64/metaprojects/combo/stable/lib:/cvmfs/icecube.opensciencegrid.org/py3-v4.1.0/RHEL_7_x86_64/metaprojects/combo/stable/cernroot/lib:/cvmfs/icecube.opensciencegrid.org/py3-v4.1.0/RHEL_7_x86_64/metaprojects/combo/stable/lib/tools:/cvmfs/icecube.opensciencegrid.org/py3-v4.1.0/RHEL_7_x86_64/spack/opt/spack/linux-centos7-x86_64/gcc-9.2.0spack/geant4-10.04-jtcfrttq4krvbvabfjtkq7itsf3dieqg/lib64
export LIBRARY_PATH=
export PYTHONPATH=/cvmfs/icecube.opensciencegrid.org/py3-v4.1.0/RHEL_7_x86_64/metaprojects/combo/stable/lib
export I3_SRC=/cvmfs/icecube.opensciencegrid.org/py3-v4.1.0/RHEL_7_x86_64/metaprojects/combo/stable
export I3_BUILD=/cvmfs/icecube.opensciencegrid.org/py3-v4.1.0/RHEL_7_x86_64/metaprojects/combo/stable
export I3_PORTS=
export I3_TESTDATA=/cvmfs/icecube.opensciencegrid.org/py3-v4.1.0/../data/i3-test-data-svn/trunk
export ROOTSYS=/cvmfs/icecube.opensciencegrid.org/py3-v4.1.0/RHEL_7_x86_64/metaprojects/combo/stable/cernroot
export LHAPDF_DATA_PATH=/cvmfs/sft.cern.ch/lcg/external/lhapdfsets/current

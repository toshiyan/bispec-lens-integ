#!/bin/sh

cwd=$(pwd)

# define functions
makeF90()
{
  echo '---- compiling F90 ' ${1} '----'
  cd F90/${1}
  make -f Makefile clean
  make -f Makefile
  make -f Makefile install
  make -f Makefile clean
  cd ${cwd}
}

makef2py()
{
  echo '----' ${1} '----'
  cd ${1}
  ./makefile.sh all
  cd ${cwd}
}

if [ ${1} = "clean" ]; then
  rm -rf wrap/*.so
  rm -rf F90/lib/*.a
  rm -rf F90/mod/*.mod
fi

# compile f90 sources
if [ ${1} = "F90" -o ${1} = "all" ]; then
  rm -rf F90/lib/*.a
  rm -rf F90/mod/*.mod
  makeF90 src
fi

# create python modules
if [ ${1} = "f2py" -o ${1} = "all" ]; then
  makef2py bispec
fi


#!/usr/bin/env sh

HOME=$(pwd)


echo "****************************************\n"
echo "*********** 3D BFM Design Tool **********\n"
echo "****************************************\n\n"
echo "Add Following in your .bashrc\n"
echo "\texport M2BFM=$HOME/"
echo "\tPATH=\$PATH:$HOME/executables/"
echo "\tPYTHONPATH=\$PYTHONPATH:$HOME/executables/\n\n"
echo "*************************************\n\n"

PYTHON=$(where python)

if [ "$PYTHON" != "/usr/bin/python" ]; then
  echo "IMPORTANT !!! \n Python directory is different than that used in the bin file. \nPlease change the header in bin/MakeBlade.py file to:\n"
  echo $PYTHON
  echo "\n"
fi

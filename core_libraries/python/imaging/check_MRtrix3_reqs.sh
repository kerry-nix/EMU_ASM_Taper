#!/bin/sh

#### Check if the required programs exist

# Python minimum requirement
pymajor=2
pyminor=6

if command -v python; then
    version=$(python --version | sed 's/^[^ ]* //')
    major=$(echo $version | cut -d '.' -f 1)
    minor=$(echo $version | cut -d '.' -f 2)
    
    if [[ "$major">"$pymajor" ]]; then
        :
    elif [[ "$major"=="$pymajor" ]] && [[ "$minor">="$pyminor" ]]; then
        :
    else
        echo "Python version $major.$minor is not supported by MRtrix3. Please upgrade Python to at least 2.6,change your default Pythonpath to a supported version, or run in the correct environment."
        exit
    fi
else
    echo "Python not found. Please install Python 2.6+, change your default Pythonpath to a supported version, or run in the correct environment."
    exit
fi
echo "Python requirements met."

# C++ Compiler minimum requirements
cppmajor=4
cppminor=9

if command -v g++; then
    version=$(g++ --version | sed 's/^[^ ]* //')
    major=$(echo $version | cut -d '.' -f 1)
    minor=$(echo $version | cut -d '.' -f 2)

    if [[ "$major">"$cppmajor" ]]; then
        :
    elif [[ "$major"=="$cppmajor" ]] && [[ "$minor">="$cppminor" ]]; then
        :
    else
        echo "g++ version $major.$minor is not supported by MRtrix3. Please upgrade g++ to at least 4.9,change your default path to a supported version, or run in the correct environment."
        exit
    fi
else
    echo "g++ not found. Please install g++ 4.9+, change your default path to a supported version, or run in the correct environment."
    exit
fi
echo "g++ requirements met."

# Eigen Compiler minimum requirements
eigenmajor=3
eigenminor=2
eigenfix=8

:<<'comment'
if command -v eigen; then
    version=$(eigen --version | sed 's/^[^ ]* //')
    major=$(echo $version | cut -d '.' -f 1)
    minor=$(echo $version | cut -d '.' -f 2)
    fix=$(echo $version | cut -d '.' -f 3)

    if [[ "$major">"$eigenmajor" ]]; then
        :
    elif [[ "$major"=="$eigenmajor" ]] && [[ "$minor">="$eigenminor" ]]; then
        :
    elif [[ "$major"=="$eigenmajor" ]] && [[ "$minor"=="$eigenminor" ]] && [[ "$fix">="$eigenfix" ]]; then
        :
    else
        echo "eigen version $major.$minor.$$fix is not supported by MRtrix3. Please upgrade eigen to at least 3.2.8,change your default path to a supported version, or run in the correct environment."
        exit
    fi
else
    echo "eigen not found. Please install eigen 3.2.8+, change your default path to a supported version, or run in the correct environment."
    exit
fi
echo "eigen requirements met."
comment

echo "Is Eigen >=$eigenmajor.$eigenminor.$eigenfix installed? Y/n"
read eigeninput
if [[ $eigeninput=='Y' ]]; then
    :
else
    echo "Please install eigen 3.2.8+, change your default path to a supported version, or run in the correct environment."
    exit
fi

# OpenGL minimum requirements
openglmajor=3
openglminor=3

echo "Is OpenGL >= $openglmajor.$openglminor installed? Y/n"
read openglinput
if [[ $openglinput=='Y' ]]; then
    :
else
    echo "Please install OpenGL >=3.3, change your default path to a supported version, or run in the correct environment."
    exit
fi

# Qt backend requirements
qtmajor=4
qtminor=8

echo "Is QT >= $qtmajor.$qtminor installed? Y/n"
read qtinput
if [[ $qtinput=='Y' ]]; then
    :
else
    echo "Please install QT >=4.8, change your default path to a supported version, or run in the correct environment."
    exit
fi

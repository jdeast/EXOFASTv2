#!/bin/bash

# this script will use git to install exofastv2 and its dependences, and test it on
# linux or mac with a zsh, tcsh, or bash shell
# it assumes IDL and git are installed and in your path
# NOTE: this will pollute your profile file (e.g., ~/.bashrc) if run multiple times


# make user's idl directory (if it doesn't already exist)
mkdir -p $HOME/idl
cd $HOME/idl

# install dependencies if they don't already exist
if command -v git &> /dev/null; then
    if [ ! -d EXOFASTv2 ]; then git clone https://github.com/jdeast/EXOFASTv2.git;   fi
    if [ ! -d IDLAstro ];  then git clone https://github.com/wlandsman/IDLAstro.git; fi
    if [ ! -d coyote ];    then git clone https://github.com/idl-coyote/coyote.git;  fi
else
    echo "git not found. Install git (or add it to $PATH) and rerun this code."
    exit 0
fi

# determine which shell we're using
#MYSHELL=`echo $0` 
if [ $SHELL == '/bin/bash' ]; then
    PROFILE_FILE=${HOME}/.bashrc
elif [ $SHELL == '/bin/zsh' ]; then
    PROFILE_FILE=${HOME}/.zshrc
elif [ $SHELL == '/bin/tcsh' ]; then
    PROFILE_FILE=${HOME}/.tcshrc
else
    echo "shell not recognized:" $SHELL
fi

# set environment variables
echo '' >> $PROFILE_FILE
echo '# added by $EXOFAST_PATH/setup.sh on' `date +%Y-%m-%dT%H:%M:%S` >> $PROFILE_FILE
if [ $SHELL == '/bin/bash' ] || [ $SHELL == '/bin/zsh' ]; then
    export 'EXOFAST_PATH="${HOME}/idl/EXOFASTv2/"' >> $PROFILE_FILE
    echo 'export EXOFAST_PATH="${HOME}/idl/EXOFASTv2/"' >> $PROFILE_FILE
    if [ -z "$IDL_PATH" ]; then
 	echo 'export IDL_PATH="<IDL_DEFAULT>:+${EXOFAST_PATH}"' >> $PROFILE_FILE
    else
 	# otherwise, append EXOFAST_PATH to your IDL_PATH
 	echo 'export IDL_PATH="${IDL_PATH}:+${EXOFAST_PATH}"' >> $PROFILE_FILE
    fi
    # add dependencies to IDL_PATH
    echo 'export IDL_PATH="${IDL_PATH}:+${HOME}/idl/IDLAstro"' >> $PROFILE_FILE
    echo 'export IDL_PATH="${IDL_PATH}:+${HOME}/idl/coyote"' >> $PROFILE_FILE
elif [ $SHELL == '/bin/tcsh' ]; then
    setenv EXOFAST_PATH "${HOME}/idl/EXOFASTv2/"
    if (! $?IDL_PATH) then
 	echo 'setenv IDL_PATH "<IDL_DEFAULT>:+{$EXOFAST_PATH}"' >> $PROFILE_FILE
    else
 	# otherwise, append EXOFAST_PATH to your IDL_PATH
 	echo 'setenv IDL_PATH "${IDL_PATH}:+${EXOFAST_PATH}"' >> $PROFILE_FILE
    fi
    # add dependencies to IDL_PATH    
    echo 'setenv IDL_PATH "${IDL_PATH}:+${HOME}/idl/IDLAstro"' >> $PROFILE_FILE
    echo 'setenv IDL_PATH "${IDL_PATH}:+${HOME}/idl/coyote"' >> $PROFILE_FILE
fi

# source the PROFILE_FILE to add environment variables to this script
source $PROFILE_FILE

# if IDL is not in the path, try to find it and add it to the path
if [ command -v idl &> /dev/null ]; then
    idl -e fithat3, nthread=1, maxsteps=100
elif [ alias idl &> /dev/null ]; then
    idl -e fithat3, nthread=1, maxsteps=100    
else
    echo "idl not found. Install idl (or add it to PATH)"
fi

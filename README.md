EXOFASTv2 -- Jason Eastman (jason.eastman@cfa.harvard.edu) An
exoplanet transit and radial velocity fitting software package in IDL
If you use this in a publication, please cite:
http://adsabs.harvard.edu/abs/2017ascl.soft10003E

A tutorial with exercises can be found here:
https://docs.google.com/document/d/1H-HMe1No5B4JE93V9kSEW91uSIQcG2eUVScf037biZw/edit

# Installation instructions #

License-free use still requires a (free) IDL installation and runs a
pre-compiled version with the virtual machine. You still must follow
the installation instructions below.

Note the IDL Astronomy library is required. If you don't already have
it, install it from here: https://github.com/wlandsman/IDLAstro

This package is best installed with git

  cd $HOME/idl
  git clone https://github.com/jdeast/EXOFASTv2.git

NOTE: EXOFASTv2 logs the version number used for each fit if installed
with git and git can be invoked via "git".

define environment variables (bash shell, e.g., .bashrc)

  EXOFAST_PATH="$HOME/idl/EXOFASTv2/" ; export EXOFAST_PATH
  # if IDL_PATH is not defined, add EXOFAST_PATH and subdirectories to the default IDL path
  if [ -z "$IDL_PATH" ]; then 
     IDL_PATH="<IDL_DEFAULT>:+${EXOFAST_PATH}" ; export IDL_PATH
  else 
     # otherwise, append EXOFAST_PATH and all subdirectories to your IDL_PATH
     IDL_PATH="${IDL_PATH}:+${EXOFAST_PATH}" ; export IDL_PATH
  fi

--- OR ---

define environment variables (c shell, e.g., .tcshrc)

  setenv EXOFAST_PATH "${HOME}/idl/EXOFASTv2/"
  # if IDL_PATH is not defined, add EXOFAST_PATH and subdirectories to the default IDL path
  if ("$IDL_PATH" == "") then 
     setenv IDL_PATH "<IDL_DEFAULT>:+${EXOFAST_PATH}"
  else
     # otherwise, append EXOFAST_PATH and all subdirectories to your IDL_PATH
     setenv IDL_PATH "${IDL_PATH}:+${EXOFAST_PATH}"
  endif


NOTE: If you have used the old version of EXOFAST, you must remove it
from your IDL_PATH to run correctly.

To test your installation, try the HAT-3b example:

  cd $EXOFAST_PATH/examples/hat3
  idl -e "fithat3"

TIP 1: make your terminal wide so the MCMC updates don't spam the screen

TIP 2: If you don't care about the results, run a very short fit by
setting maxsteps=100 (idl -e "fithat3, maxsteps=100"). This will give
you a very imprecise/unreliable answer, but allow you to check your
installation in ~5 seconds instead of ~10 minutes.

It will generate many output files
($EXOFAST_PATH/examples/hat3/HAT-P-3b.Torres.*) (see
$EXOFAST_PATH/exofastv2.pro for an explanation of outputs). The last
file it generates is HAT-3b.Torres.chain.ps. If that is generated
without error, you're good to go!

To get future updates, simply type

  cd $EXOFAST_PATH
  git pull

# Troubleshooting #

I try hard to test thoroughly before pushing new code (but I'm not
perfect!). If it does not compile or you get a syntax error, it is
very likely a problem with your setup. The most likely reasons are:

1) Your IDL_PATH or EXOFAST_PATH environment variables are not set up
correctly. From a terminal, type "echo $IDL_PATH" (it should include
EXOFASTv2) and "echo $EXOFAST_PATH" (it should point to your
installation) to check.

2) You have missing dependencies (e.g., IDL astronomy library, coyote
library)

3) You have programs with the same name with a higher precedence in
your IDL path. Renaming or moving your version will fix it, but please
send me an email if the conflicting code is a library routine. I will
rename the EXOFASTv2 version to avoid conflicts with others.

4) You have an incompatible version of IDL. EXOFASTv2 has been built
and tested on linux with IDL 8.5. I also typically perform a cursory
test with a linux installation of IDL 6.4 (circe 2007). I am not aware
of any incompatibility for any platform (Windows, Mac, Linux) or IDL
versions newer than 5.0, but it has not been tested on anything older
than 6.4. If you find any incompatibilites on any version, IDL 5.0
(circe 1997) or later, please let me know. This code relies heavily on
pointers and structures, which were introduced in IDL 5.0. Older
versions will never be supported.

Note 1: The latest IDL version can be installed for free and EXOFASTv2
can be run within a virtual machine without a license.

Note 2: The HAT-3 example runs to completion in GDL but many things are
sub-optimal:

   a) The results are unverified and many features are untested 
   b) Limb darkening models outside of the grid are rejected a priori.
   c) Multi-page postscript files (chains, pdfs, models) are not supported in GDL.
   d) The covariance plot is disabled
   e) It is ~3.5x slower
   f) Error messages spam the screen

If you're interested in making it work better with GDL, please contact
me.

5) I have introduced a bug. Even if your problem is not a bug, if
you've read the documentation, given it some thought, and still can't
figure it out, it's probably at least a failure in documentation. Send
me an email.

# Tips, Warnings and Caveats #

Other examples are available for various use cases, which are intended
to be templates for various types of fits, and includes an example of
running EXOFASTv2 without an IDL license. See
$EXOFAST_PATH/examples/README for more information.

Error checking is not thorough. You may encounter cryptic error
messages and strange failure modes if you use it in a way it wasn't
intended. Do not stray far from the examples blindly. 

Do not ignore the warnings about the Gelman Rubin statistic without
thoroughly inspecting the PDFs and chains.

If you're stuck, feel free to ask me for guidance.

This is a BETA version. EXPECT BUGS!!! And please report any
unexpected behavior.

It is not fully documented. Please don't hesitate to email me with
questions. See the $EXOFAST_PATH/examples directory for templates to
get started on your own fits.

Major releases or bug fixes will be announced on twitter
(@exofastupdates)



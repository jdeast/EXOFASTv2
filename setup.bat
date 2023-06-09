@setx MODELDIR %HOME%\modeling\
@setx EXOFAST_PATH %HOME%\idl\EXOFASTv2
@setx IDL_PATH "<IDL_DEFAULT>:+$HOME/idl"

:: This is hardcoded to 64 bit IDL v8.8. Change to your actual path and uncomment before executing
@setx PATH "%PROGRAMFILES%\Harris\IDL88\bin\bin.x86_64\"

@cd %HOME%\idl
if exist EXOFASTv2\ (
  echo "EXOFASTv2 already installed"
) else (
  @git clone https://github.com/jdeast/EXOFASTv2.git
)

if exist IDLAstro\ (
  echo "IDLAstro already installed"
) else (
  @git clone https://github.com/wlandsman/IDLAstro.git
)

if exist coyote\ (
  echo "coyote already installed"
) else (
  @git clone https://github.com/idl-coyote/coyote.git
)

@%EXOFAST_PATH%\examples\hat3_nolicense\fithat3.bat

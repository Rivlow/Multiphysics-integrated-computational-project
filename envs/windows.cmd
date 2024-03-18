@echo off
:: This file opens a terminal which allows you to compile the code with
:: a 64-bit mingw compiler
::     https://mingw-w64.org/
::     https://www.msys2.org/
::
:: HOW TO USE THIS FILE?
::   [check the PATHs below]
::   [run this file]
::   mkdir build
::   cd build
::   cmake ..
::   make
::   ctest
::
:: How to clean the "build" folder using cmd line?
::   cd build
::   rd /q /s .

echo setting math0471 environment...

:: set the location of mingw compiler / zlib / cmake
set COMPILERPATH=%USERPROFILE%\mingw64\bin
set ZLIB=%~dp0..\lib\zlib

:: use ninja (UTF-8 chars in USERPROFILE => OK)
set MAKESYSTEM=%~dp0ninja
:: use mingw32-make (UTF-8 chars in USERPROFILE => KO)
:: (uncomment the following line if you want to use classical Makefiles)
:: set MAKESYSTEM=%~dp0make

:: look for utilities

FOR /F "tokens=* USEBACKQ" %%F IN (`where powershell`) DO ( SET POWERSHELLPATH=%%F )
ECHO POWERSHELLPATH   = %POWERSHELLPATH%
FOR /F "tokens=* USEBACKQ" %%F IN (`where bitsadmin`) DO ( SET BITSADMINPATH=%%F )
ECHO BITSADMINPATH    = %BITSADMINPATH%
FOR /F "tokens=* USEBACKQ" %%F IN (`where cmake`) DO ( SET CMAKEPATH=%%F )
ECHO CMAKEPATH        = %CMAKEPATH%
FOR /F "tokens=* USEBACKQ" %%F IN (`where git`) DO ( SET GITPATH=%%F )
ECHO GITPATH          = %GITPATH%

FOR %%I IN ("%POWERSHELLPATH%") DO SET "POWERSHELLFOLDER=%%~dpI"
ECHO POWERSHELLFOLDER = %POWERSHELLFOLDER%
FOR %%I IN ("%BITSADMINPATH%") DO SET "BITSADMINFOLDER=%%~dpI"
ECHO BITSADMINFOLDER  = %BITSADMINFOLDER%
FOR %%I IN ("%CMAKEPATH%") DO SET "CMAKEFOLDER=%%~dpI"
ECHO CMAKEFOLDER      = %CMAKEFOLDER%
FOR %%I IN ("%GITPATH%") DO SET "GITFOLDER=%%~dpI"
ECHO GITFOLDER        = %GITFOLDER%

:: clears the PATH!
:: avoids things to be found in paths such as "c:\Strawberry\perl\bin", 
:: or "c:\local"!

set "PATH=%COMPILERPATH%"
set "PATH=%PATH%;%CMAKEFOLDER%"
set "PATH=%PATH%;%BITSADMINFOLDER%"
set "PATH=%PATH%;%POWERSHELLFOLDER%"
set "PATH=%PATH%;%GITFOLDER%"

:: add current folder to PATH for cmake/make aliases
set PATH=%MAKESYSTEM%;%PATH%
:: where is zlib.dll ?
set PATH=%ZLIB%\bin;%PATH%
:: where is zlib.h ?
set INCLUDE=%ZLIB%\include
:: where is zlib.lib ?
set LIB=%ZLIB%\lib

:: clear local vars
set COMPILERPATH=
set ZLIB=
set MAKESYSTEM=
set POWERSHELLPATH=
set BITSADMINPATH=
set CMAKEPATH=
set GITPATH=
set POWERSHELLFOLDER=
set BITSADMINFOLDER=
set CMAKEFOLDER=
set GITFOLDER=

:: status
echo PATH    = %PATH%
echo INCLUDE = %INCLUDE%
echo LIB     = %LIB%

:: open terminal
CD /d "%~dp0"
CD ..
CD src/build
%comspec% /K



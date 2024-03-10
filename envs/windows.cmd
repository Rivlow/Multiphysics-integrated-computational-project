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

echo setting MinGW64 environment...

:: set the location of mingw compiler / zlib
set COMPILERPATH=%USERPROFILE%\mingw64\bin
set ZLIB=%~dp0..\lib\zlib
:: use ninja (UTF-8 chars in USERPROFILE => OK)
set MAKESYSTEM=%~dp0ninja
:: use mingw32-make (UTF-8 chars in USERPROFILE => KO)
:: (uncomment the following line if you want to use classical Makefiles)
:: set MAKESYSTEM=%~dp0make

:: Look for MinGW.
IF NOT EXIST "%COMPILERPATH%" (
    ECHO   - compiler NOT found in %COMPILERPATH%!!
    @REM PAUSE
    @REM EXIT /B
) ELSE (
    ECHO   - compiler found in %COMPILERPATH%.
    set "PATH=%COMPILERPATH%;%PATH%"
)

IF NOT EXIST "%ZLIB%" (
    ECHO   - zlib NOT found in %ZLIB%!!
    @REM PAUSE
    @REM EXIT /B
) ELSE (
    ECHO   - zlib found in %ZLIB%.
)

:: add current folder to PATH for cmake/make aliases
set PATH=%MAKESYSTEM%;%PATH%
:: where is zlib.dll ?
set PATH=%ZLIB%\bin;%PATH%
:: where is zlib.h ?
set INCLUDE=%ZLIB%\include;%INCLUDE%
:: where is zlib.lib ?
set LIB=%ZLIB%\lib;%LIB%

:: clear local vars
set COMPILERPATH=
set ZLIB=
set MAKESYSTEM=

:: open terminal
CD /d "%~dp0"
CD ..
CD src/build
%comspec% /K



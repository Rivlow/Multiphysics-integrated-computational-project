::@echo off
:: download zlib & extract & build
:: https://eigen.tuxfamily.org/

setlocal
set version=1.3.1
set folder=zlib-%version%
set file=%folder%.zip

:: download sources from github
del %file%  >nul 2>&1
bitsadmin /transfer get_zlib /dynamic /download /priority foreground "https://github.com/madler/zlib/archive/refs/tags/v%version%.zip" "%CD%\%file%"

:: unzip sources
rd /Q/S %folder%  >nul 2>&1
powershell.exe -nologo -noprofile -command "& { Add-Type -A 'System.IO.Compression.FileSystem'; [IO.Compression.ZipFile]::ExtractToDirectory('%CD%\%file%', '%CD%'); }"
del %file%  >nul 2>&1

:: build from source
cd %folder%
mkdir build
cd build

cmake.exe -G"Ninja" -DCMAKE_INSTALL_PREFIX=..\..\zlib -DCMAKE_BUILD_TYPE=Release ..
cmake.exe --build . --target install
cd ..\..

:: remove src folder
rd /Q/S %folder% >nul 2>&1

:: clear local vars
endlocal

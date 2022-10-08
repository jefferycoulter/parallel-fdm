@echo OFF

rem get current directory
set CURRENT_DIR=%cd%
pushd %CURRENT_DIR%
rem run cmake
cmake -S %CURRENT_DIR% -B %CURRENT_DIR%\build ^
        -G "Unix Makefiles" ^
        -DCMAKE_C_COMPILER=clang ^
        -DCMAKE_MAKE_PROGRAM="C:\\Program Files (x86)\\GnuWin32\\bin\\make.exe" ^
        -DCMAKE_BUILD_TYPE=Debug
cd build
rem tell make to run on multiple threads if possible
if "%MAKEFLAGS%" == "" set MAKEFLAGS=-j%NUMBER_OF_PROCESSORS%
rem run make
make
popd

pause
rem Build GLPK DLL with Microsoft Visual Studio Community 2015

rem NOTE: Make sure that HOME variable specifies correct path
set VCVARSALL_PATH="C:\Program Files (x86)\Microsoft Visual Studio\2017\Community\VC\Auxiliary\Build"
set NMAKE_PATH="C:\Program Files (x86)\Microsoft Visual Studio\2017\Community\VC\Tools\MSVC\14.16.27023\bin\Hostx64\x64"

call %VCVARSALL_PATH%\vcvarsall.bat x64
copy config_VC config.h
%NMAKE_PATH%\nmake.exe /f Makefile_VC
%NMAKE_PATH%\nmake.exe /f Makefile_VC check

pause
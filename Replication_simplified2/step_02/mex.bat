@IF "%VSINSTALLDIR%"=="" call "C:\Program Files (x86)\Microsoft Visual Studio 11.0\VC\bin\x86_amd64\vcvarsx86_amd64.bat"

@set MATLAB=C:\Program Files\MATLAB\R2013b

cl /I "%MATLAB%\extern\include" /LD /EHsc /O2 /Fe%~n1.mexw64 %1 libmx.lib libmex.lib libmwblas.lib /link "/EXPORT:mexFunction" "/LIBPATH:%MATLAB%\extern\lib\win64\microsoft"

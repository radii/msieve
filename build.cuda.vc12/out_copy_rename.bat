@echo off
if not exist %1 goto nofile
if exist %2 goto next

echo creating directory %2
md %2 > nul

:next
rem strip quotes if present
set str=%2
for /f "useback tokens=*" %%a in ('%str%') do set str=%%~a

rem add a backslash if the output directory lacks one
set str=%str:~-1%
if "%str%" == "\" (set outf=%2%3) else (set outf=%2\%3)

echo copying %1 to %outf% (if not present or changed)
if not exist "%outf%" goto copy

rem don't overwrite if output exists and is not changed
fc %1 %outf% > nul && if not %errorlevel 1 goto exit
echo overwriting %outf% with %1

:copy
copy %1 %outf% > nul
goto exit

:nofile
echo %1 not found

:exit

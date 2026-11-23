@echo off
setlocal
cd /d "%~dp0"

rem Double-click friendly. With no arguments the console opens on the default
rem region and resumes from the last checkpointed block; arguments override the
rem defaults one by one. Nothing is scanned until Start is pressed on the page.
set DLO=1.4815e18
set DHI=2e18
set THREADS=6
set OUT=runs\k114

if /i "%~1"=="help" goto usage
if /i "%~1"=="-h" goto usage
if /i "%~1"=="--help" goto usage

if not "%~1"=="" set DLO=%~1
if not "%~2"=="" set DHI=%~2
if not "%~3"=="" set THREADS=%~3
if not "%~4"=="" set OUT=%~4

if not exist "target\release\k114.exe" (
  echo engine not built - running cargo build --release ...
  cargo build --release
  if errorlevel 1 exit /b 1
)

python app\server.py --dlo %DLO% --dhi %DHI% --threads %THREADS% --out %OUT% --open
if errorlevel 1 (
  echo.
  echo console exited with an error - see the messages above.
  pause
)
exit /b 0

:usage
echo usage: console.bat [dlo] [dhi] [threads] [out-dir]
echo.
echo   no arguments  scan d in [1.4815e18, 2e18), 6 threads, out runs\k114,
echo                 resuming after the last recorded block
echo   example       console.bat 1e12 1.00000016777216e12 4 runs\scratch
echo.
echo R = ^|z^|/d is fixed at [3.4, 64]. Press Start on the page to scan and Stop
echo to checkpoint. The upper edge 2e18 is a formality: at the measured
echo throughput it is not reachable, so the region never reports "complete".
exit /b 2

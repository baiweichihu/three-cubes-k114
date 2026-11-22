@echo off
setlocal
cd /d "%~dp0"

if "%~1"=="" (
  echo usage: console.bat ^<dlo^> ^<dhi^> [threads] [out-dir]
  echo   example: console.bat 1.4815e18 1.48152e18 2 runs\region-a
  echo.
  echo The console scans d in [dlo, dhi) with R = ^|z^|/d in [3.4, 64],
  echo keeps 2 worker threads by default, and resumes from the last
  echo checkpointed block on every Start.
  exit /b 2
)

set DLO=%~1
set DHI=%~2
set THREADS=%~3
if "%THREADS%"=="" set THREADS=2
set OUT=%~4
if "%OUT%"=="" set OUT=runs\k114

if not exist "target\release\k114.exe" (
  echo engine not built - running cargo build --release ...
  cargo build --release
  if errorlevel 1 exit /b 1
)

python app\server.py --dlo %DLO% --dhi %DHI% --threads %THREADS% --out %OUT% --open

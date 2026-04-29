# Copy the gfortran runtime DLLs from your gfortran installation into
# both lib/ (as a backup) and bin/ (so interfaz_remix.exe can be run /
# distributed without a gfortran installation on the target machine).
#
# Usage: ./copy_dlls.ps1

$ErrorActionPreference = 'Stop'

Write-Host "Copying gfortran runtime DLLs..." -ForegroundColor Green

# Locate gfortran on PATH
$gfortranCmd = Get-Command gfortran -ErrorAction SilentlyContinue
if (-not $gfortranCmd) {
    Write-Error "gfortran not found in PATH. Install MinGW/MSYS2/TDM-GCC and add it to PATH, or use copy_local_dlls.ps1 instead."
    exit 1
}

$gccPath = Split-Path $gfortranCmd.Source
Write-Host "Found gfortran at: $gccPath"

# Resolve repository paths relative to this script
$root = Split-Path -Parent $MyInvocation.MyCommand.Path
$libDir = Join-Path $root 'lib'
$binDir = Join-Path $root 'bin'

foreach ($d in @($libDir, $binDir)) {
    if (-not (Test-Path $d)) {
        New-Item -ItemType Directory -Force -Path $d | Out-Null
    }
}

# Required runtime DLLs
$dlls = @(
    'libgfortran_64-5.dll',
    'libgcc_s_seh_64-1.dll',
    'libwinpthread_64-1.dll',
    'libquadmath_64-0.dll'
)

$copied  = 0
$missing = @()
foreach ($dll in $dlls) {
    $src = Join-Path $gccPath $dll
    if (Test-Path $src) {
        Copy-Item $src (Join-Path $libDir $dll) -Force
        Copy-Item $src (Join-Path $binDir $dll) -Force
        Write-Host "Copied $dll -> lib/ and bin/" -ForegroundColor Green
        $copied++
    } else {
        Write-Warning "Not found in $gccPath : $dll"
        $missing += $dll
    }
}

Write-Host ""
Write-Host "Copied $copied / $($dlls.Count) DLLs."
if ($missing.Count -gt 0) {
    Write-Host "Missing: $($missing -join ', ')" -ForegroundColor Yellow
    exit 1
}

Write-Host "interfaz_remix.exe is now portable!" -ForegroundColor Green
exit 0

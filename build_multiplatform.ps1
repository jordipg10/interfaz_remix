# Build the interfaz_remix executable for multiple platforms.
#
# - Windows: copies the existing bin/interfaz_remix.exe (and its runtime
#   DLLs) into dist/.
# - Linux:   builds via Docker if Dockerfile.linux exists and Docker is
#   available; otherwise prints instructions.
#
# Usage: ./build_multiplatform.ps1

$ErrorActionPreference = 'Stop'

Write-Host "Building interfaz_remix for multiple platforms..." -ForegroundColor Green

# Resolve repo root and key folders
$root    = Split-Path -Parent $MyInvocation.MyCommand.Path
$binDir  = Join-Path $root 'bin'
$distDir = Join-Path $root 'dist'

# Output directory
if (-not (Test-Path $distDir)) {
    New-Item -ItemType Directory -Path $distDir | Out-Null
}

# ---------------------------------------------------------------------------
# Linux build (Docker)
# ---------------------------------------------------------------------------
$dockerfile     = Join-Path $root 'Dockerfile.linux'
$dockerAvailable = Get-Command docker -ErrorAction SilentlyContinue

if ($dockerAvailable -and (Test-Path $dockerfile)) {
    Write-Host "Building Linux executable using Docker..." -ForegroundColor Yellow
    docker build -f $dockerfile -t interfaz-remix-linux $root
    docker run --rm -v "${distDir}:/output" interfaz-remix-linux
    Write-Host "Linux executable built." -ForegroundColor Green
} else {
    if (-not $dockerAvailable) {
        Write-Warning "Docker not found. Skipping Linux build."
    }
    if (-not (Test-Path $dockerfile)) {
        Write-Warning "Dockerfile.linux not found. Skipping Linux build."
    }
    Write-Host "To build for Linux:" -ForegroundColor Yellow
    Write-Host "  1. Install Docker Desktop and add a Dockerfile.linux at the repo root, OR" -ForegroundColor Yellow
    Write-Host "  2. Compile natively on a Linux system (see BUILD_GUIDE.md)." -ForegroundColor Yellow
}

# ---------------------------------------------------------------------------
# Windows build (uses the already-linked binary in bin/)
# ---------------------------------------------------------------------------
Write-Host "Packaging Windows executable..." -ForegroundColor Yellow

$winExe = Join-Path $binDir 'interfaz_remix.exe'
if (Test-Path $winExe) {
    Copy-Item $winExe (Join-Path $distDir 'interfaz_remix.exe') -Force
    Write-Host "Copied interfaz_remix.exe -> dist/" -ForegroundColor Green
} else {
    Write-Warning "Windows executable not found at $winExe. Run the 'rebuild' VS Code task first (see BUILD_GUIDE.md)."
}

# Copy DLLs for Windows portability
$dlls = @(
    'libgfortran_64-5.dll',
    'libgcc_s_seh_64-1.dll',
    'libwinpthread_64-1.dll',
    'libquadmath_64-0.dll'
)

$copiedDlls = 0
foreach ($dll in $dlls) {
    $srcDll = Join-Path $binDir $dll
    if (Test-Path $srcDll) {
        Copy-Item $srcDll $distDir -Force
        $copiedDlls++
    } else {
        Write-Warning "DLL not present in bin/: $dll (run copy_local_dlls.ps1 or copy_dlls.ps1)"
    }
}

if ($copiedDlls -gt 0) {
    Write-Host "Copied $copiedDlls DLL files." -ForegroundColor Green
}

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
Write-Host ""
Write-Host "Build complete. Contents of $distDir :" -ForegroundColor Green
Get-ChildItem $distDir | ForEach-Object {
    Write-Host "  $($_.Name)" -ForegroundColor Cyan
}

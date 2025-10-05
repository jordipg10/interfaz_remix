# Build script for cross-platform executables
# Usage: .\build_multiplatform.ps1

Write-Host "Building interfaz for multiple platforms..." -ForegroundColor Green

# Create output directory
if (-not (Test-Path "dist")) {
    New-Item -ItemType Directory -Path "dist"
}

# Build for Linux using Docker
Write-Host "Building Linux executable..." -ForegroundColor Yellow
docker build -f Dockerfile.linux -t interfaz-linux .
docker run --rm -v "${PWD}\dist:/output" interfaz-linux

# Build for Windows (current platform)
# Build script for cross-platform executables
# Usage: .\build_multiplatform.ps1

Write-Host "Building interfaz for multiple platforms..." -ForegroundColor Green

# Create output directory
if (-not (Test-Path "dist")) {
    New-Item -ItemType Directory -Path "dist"
}

# Check if Docker is available
$dockerAvailable = Get-Command docker -ErrorAction SilentlyContinue

if ($dockerAvailable) {
    # Build for Linux using Docker
    Write-Host "Building Linux executable using Docker..." -ForegroundColor Yellow
    docker build -f Dockerfile.linux -t interfaz-linux .
    docker run --rm -v "${PWD}\dist:/output" interfaz-linux
    Write-Host "Linux executable built" -ForegroundColor Green
} else {
    Write-Warning "Docker not found. Skipping Linux build."
    Write-Host "To build for Linux:" -ForegroundColor Yellow
    Write-Host "  1. Install Docker Desktop" -ForegroundColor Yellow
    Write-Host "  2. Re-run this script" -ForegroundColor Yellow
    Write-Host "  OR compile natively on a Linux system" -ForegroundColor Yellow
}

# Build for Windows (current platform)
Write-Host "Building Windows executable..." -ForegroundColor Yellow
if (Test-Path "exe\interfaz.exe") {
    Copy-Item "exe\interfaz.exe" "dist\interfaz_windows.exe"
    Write-Host "Windows executable copied" -ForegroundColor Green
} else {
    Write-Warning "Windows executable not found. Run build task first."
}

# Copy DLLs for Windows
$dlls = @("libgfortran_64-5.dll", "libgcc_s_seh_64-1.dll", "libwinpthread_64-1.dll", "libquadmath_64-0.dll")
$copiedDlls = 0
foreach ($dll in $dlls) {
    if (Test-Path "exe\$dll") {
        Copy-Item "exe\$dll" "dist\"
        $copiedDlls++
    }
}

if ($copiedDlls -gt 0) {
    Write-Host "Copied $copiedDlls DLL files" -ForegroundColor Green
}

Write-Host "Build complete! Check the dist folder for executables." -ForegroundColor Green
if (Test-Path "dist") {
    Get-ChildItem "dist" | ForEach-Object {
        Write-Host "  $($_.Name)" -ForegroundColor Cyan
    }
}
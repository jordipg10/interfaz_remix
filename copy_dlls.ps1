# Simple DLL copy script for interfaz.exe portability

Write-Host "Copying GFortran runtime DLLs..." -ForegroundColor Green

# Find gfortran path
$gfortranCmd = Get-Command gfortran -ErrorAction SilentlyContinue
if (-not $gfortranCmd) {
    Write-Error "gfortran not found in PATH"
    exit 1
}

$gccPath = Split-Path $gfortranCmd.Source
Write-Host "Found gfortran at: $gccPath"

# DLLs to copy
$dlls = @(
    "libgfortran_64-5.dll",
    "libgcc_s_seh_64-1.dll", 
    "libwinpthread_64-1.dll",
    "libquadmath_64-0.dll"
)

# Copy each DLL to lib folder first
$copied = 0
foreach ($dll in $dlls) {
    $src = Join-Path $gccPath $dll
    $dstLib = Join-Path "lib" $dll
    $dstExe = Join-Path "exe" $dll
    
    if (Test-Path $src) {
        # Copy to lib folder (for organization)
        Copy-Item $src $dstLib -Force
        # Copy to exe folder (for runtime)
        Copy-Item $src $dstExe -Force
        Write-Host "Copied $dll to lib/ and exe/" -ForegroundColor Green
        $copied++
    } else {
        Write-Warning "Not found: $dll"
    }
}

Write-Host "Copied $copied DLLs to lib and exe folders"
Write-Host "interfaz.exe is now portable!" -ForegroundColor Green
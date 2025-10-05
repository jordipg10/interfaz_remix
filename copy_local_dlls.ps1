Param()

$ErrorActionPreference = 'Stop'

# Resolve paths relative to this script
$root = Split-Path -Parent $MyInvocation.MyCommand.Path
$src  = Join-Path $root 'lib'
$dst  = Join-Path $root 'exe'

if (-not (Test-Path $src)) {
    Write-Error "Source folder not found: $src"
    exit 1
}

if (-not (Test-Path $dst)) {
    New-Item -ItemType Directory -Force -Path $dst | Out-Null
}

$dlls = @(
    'libgfortran_64-5.dll',
    'libquadmath_64-0.dll',
    'libgcc_s_seh_64-1.dll',
    'libwinpthread_64-1.dll'
)

$missing = @()
foreach ($d in $dlls) {
    $s = Join-Path $src $d
    if (Test-Path $s) {
        Copy-Item $s $dst -Force
        Write-Host "Copied $d" -ForegroundColor Green
    } else {
        Write-Warning "Missing $d in lib/"
        $missing += $d
    }
}

if ($missing.Count -gt 0) {
    Write-Host "Completed with warnings. Missing: $($missing -join ', ')" -ForegroundColor Yellow
} else {
    Write-Host 'All DLLs copied successfully.' -ForegroundColor Green
}

exit 0

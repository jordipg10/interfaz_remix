# Cross-Platform Build Guide

This guide explains how to build executables for different operating systems.

## Prerequisites

- General
	- A recent gfortran toolchain installed and available on PATH
	- VS Code if you want to use the provided tasks
	- Write access to create folders: `obj`, `mod`, `exe`, `dist`
- Windows
	- MinGW/TDM-GCC gfortran (64-bit)
	- Runtime DLLs are provided in this repo under `lib/`:
		- `libgfortran_64-5.dll`
		- `libquadmath_64-0.dll`
		- `libgcc_s_seh_64-1.dll`
		- `libwinpthread_64-1.dll`
- Linux
	- Install `gfortran` via your distro’s package manager
- macOS
	- Install `gcc` (gfortran) via Homebrew: `brew install gcc`

## Windows (Current Platform)

### Option 1: Use VS Code tasks
Preferred (one-shot):
- Run `build-with-local-dlls` (links dynamically and copies the four gfortran runtime DLLs from `lib/` into `exe/`).

Step-by-step:
1. Run `compile`
2. Run `link` (dynamic linking)
3. Run `copy-local-dlls` (copies DLLs from your repo `lib/` to `exe/`)

Alternative: `build-with-dlls` (uses `copy_dlls.ps1` to locate DLLs from your gfortran installation instead of the repo `lib/` folder). Use this if your `lib/` doesn’t have the DLLs.

### Option 2: Manual build
```powershell
# Compile (recommended via VS Code task due to ordering and command length):
#   Tasks: Run Task -> compile
# If compiling manually, ensure module order and length limits are handled.

# Link (dynamic; recommended)
cd obj
gfortran -ffpe-summary=none -o ../exe/interfaz.exe *.o

# Copy runtime DLLs from repo lib/ to exe/
Copy-Item ..\lib\libgfortran_64-5.dll ..\lib\libquadmath_64-0.dll ..\lib\libgcc_s_seh_64-1.dll ..\lib\libwinpthread_64-1.dll ..\exe\ -Force
```

## Linux

### Option 1: Native compilation (on Linux system)
```bash
# Install gfortran
sudo apt-get install gfortran          # Ubuntu/Debian
sudo yum install gcc-gfortran          # RHEL/CentOS  
sudo dnf install gcc-gfortran          # Fedora

# Create directories
mkdir -p obj mod exe

# Compile
cd obj
gfortran -g -c -O0 -ffpe-summary=none -ffree-line-length-none -fno-range-check -J ../mod ../src_REMIX/*.f90 ../src_interfaces/*.f90

# Link (dynamic)
gfortran -ffpe-summary=none -o ../exe/interfaz_linux *.o

# Optional: static linking (may require extra system static libs)
# gfortran -static -ffpe-summary=none -o ../exe/interfaz_linux *.o
```

### Option 2: Using Docker (from Windows)
- A `Dockerfile.linux` is required (example below). After adding it, you can run:
	```powershell
	.\build_multiplatform.ps1
	```

Example `Dockerfile.linux` (save at repo root if you want to use Docker):
```Dockerfile
FROM debian:stable-slim
RUN apt-get update && apt-get install -y gfortran && rm -rf /var/lib/apt/lists/*
WORKDIR /src
COPY . /src
RUN mkdir -p obj mod exe \
		&& cd obj \
		&& gfortran -g -c -O0 -ffpe-summary=none -ffree-line-length-none -fno-range-check -J ../mod \
				 ../src_REMIX/*.f90 ../src_interfaces/*.f90 \
		&& gfortran -ffpe-summary=none -o ../exe/interfaz_linux *.o
CMD ["/bin/bash", "-lc", "cp -v exe/interfaz_linux /output/"]
```

### Option 3: Cross-compilation (advanced)
```powershell
# Install cross-compiler
# Use build-linux VS Code task
```

## macOS

### Option 1: Native compilation (on Mac)
```bash
# Install gfortran via Homebrew
brew install gcc

# Create directories  
mkdir -p obj mod exe

# Compile
cd obj
gfortran -g -c -O0 -ffpe-summary=none -ffree-line-length-none -fno-range-check -J ../mod ../src_REMIX/*.f90 ../src_interfaces/*.f90

# Link (dynamic)
gfortran -ffpe-summary=none -o ../exe/interfaz_macos *.o

# Note: fully static linking is generally not supported on macOS.
```

### Option 2: GitHub Actions (optional)
You can add a workflow later to automate builds. No workflow is included in this repo yet.

## Distribution

### Windows
- Include `exe/interfaz.exe`
- Include these runtime DLLs next to the exe (already in your repo `lib/`):
	- `libgfortran_64-5.dll`
	- `libquadmath_64-0.dll`
	- `libgcc_s_seh_64-1.dll`
	- `libwinpthread_64-1.dll`
- You can run the task `copy-local-dlls` to place these DLLs into `exe/`, then zip the `exe/` folder for distribution. Include `DB/`, `documentation/`, and `examples/` as needed.
- A static linking variant would avoid DLLs but is not the default.

### Linux/macOS
- Single executable file when dynamically linked; ensure any needed shared libs exist on the target system
- Make executable: `chmod +x interfaz_linux`

## Available Scripts

- `build_multiplatform.ps1` - Build all platforms (requires Docker)
- `copy_dlls.ps1` - Copy runtime DLLs for portability
- `copy_local_dlls.ps1` - Copy runtime DLLs from repo `lib/` to `exe/`

## VS Code Tasks

- `compile` - Compile source files
- `link` - Link (dynamic)
- `copy-local-dlls` - Copy runtime DLLs from repo `lib/` to `exe/`
- `build-with-local-dlls` - Link + copy local DLLs (one-shot)
- `copy-dlls` - Copy runtime DLLs from your gfortran installation (alternative)
- `build-with-dlls` - Link + copy installation DLLs (alternative)
- `build-linux` - Cross-compile for Linux (requires setup)

## Troubleshooting

- VS Code task errors
	- Ensure `gfortran` is on PATH: `gfortran --version`
	- Run `clean` then try `build-with-local-dlls` again
- Missing DLL at runtime (Windows)
	- Re-run `copy-local-dlls` or ensure the four DLLs are beside `exe/interfaz.exe`
- Command length or order issues when compiling manually
	- Prefer the `compile` task (it has the correct file order)
	- If compiling outside VS Code, consider batching or a makefile
# Portable Executable Setup for interfaz.exe

This document explains how to make the `interfaz.exe` executable work on other Windows computers that don't have the GFortran compiler installed.

## The Problem

The `interfaz.exe` executable is compiled with GFortran and depends on several runtime DLL files:
- `libgfortran_64-5.dll` - Fortran runtime library
- `libgcc_s_seh_64-1.dll` - GCC runtime library  
- `libwinpthread_64-1.dll` - Threading library
- `libquadmath_64-0.dll` - Extended precision math library

These DLLs are usually installed with the compiler but are not present on computers without GFortran.

## Solutions

### Option 1: Copy Runtime DLLs (Recommended)

This is the simplest approach and has already been implemented.

#### Automatic Setup
Use one of these:
- VS Code task: `copy-local-dlls` (copies the four DLLs from your repo `lib/` into `exe/`)
- PowerShell script: `./copy_local_dlls.ps1`
- Alternative (from compiler install): `./copy_dlls.ps1` (searches your gfortran installation for the DLLs)

#### Manual Setup
1. Preferred: Copy from the repo `lib/` to the `exe/` folder (already included in this project).
   - If `lib/` is missing the files, locate your GFortran installation (e.g., `C:\TDM-GCC-64\bin\`) and copy them into `lib/` and `exe/`.
2. Copy these DLL files to both the `lib\` folder (for organization) and `exe\` folder (for runtime):
   - `libgfortran_64-5.dll`
   - `libgcc_s_seh_64-1.dll` 
   - `libwinpthread_64-1.dll`
   - `libquadmath_64-0.dll`

**Note**: The DLLs must be in the same directory as the executable (`exe\`) for the program to run properly.

#### Using VS Code Tasks
- Run the `build-with-dlls` task to build and copy DLLs in one step

### Option 2: Static Linking (optional)

You can try to embed the Fortran/GCC runtimes to avoid shipping DLLs, but this is not enabled by default.

Suggested linker flags (Windows/MinGW):
- `-static-libgfortran -static-libgcc -static-libquadmath`
   - Full `-static` is often problematic on Windows and not recommended.

How to use:
- Edit the `link` task to add the static runtime flags, or create a separate `link-static` task.
- Note: Some combinations of libraries may still require DLLs.

### Option 3: Installer Package

Create an installer that includes the executable and required DLLs:
- Use tools like NSIS, InnoSetup, or WiX
- Include all dependencies in the installer
- Handles installation and uninstallation

## Current Status

The project now organizes DLLs as follows:

**lib/ folder** (for organization):
- `libgfortran_64-5.dll` - Fortran runtime
- `libgcc_s_seh_64-1.dll` - GCC runtime
- `libwinpthread_64-1.dll` - Threading support
- `libquadmath_64-0.dll` - Extended precision math

**exe/ folder** (for runtime):
- `interfaz.exe` - The main executable
- `libgfortran_64-5.dll` - Fortran runtime (required at runtime)
- `libgcc_s_seh_64-1.dll` - GCC runtime (required at runtime)
- `libwinpthread_64-1.dll` - Threading support (required at runtime)
- `libquadmath_64-0.dll` - Extended precision math (required at runtime)

## Distribution

To distribute the executable:

1. **Zip the exe folder**: Create a zip file containing the entire `exe\` folder
2. **Test on clean system**: Verify it works on a computer without GFortran
3. **Include instructions**: Provide README with any additional setup steps

## VS Code Tasks Available

- `compile`: Compile object files into `obj/`
- `link`: Link (dynamic) into `exe/interfaz.exe`
- `copy-local-dlls`: Copy four runtime DLLs from repo `lib/` to `exe/`
- `build-with-local-dlls`: Compile → link → copy-local-dlls (one step)
- `copy-dlls`: Copy runtime DLLs from your gfortran installation (alternative)
- `build-with-dlls`: Compile + Link + copy-dlls (alternative)

## Troubleshooting

If the executable still doesn't work on other computers:

1. **Check Windows version compatibility**: Ensure target systems are compatible
2. **Verify all DLLs**: Use `dumpbin /dependents interfaz.exe` to check dependencies
3. **Missing Visual C++ runtime**: Some systems may need Visual C++ Redistributable
4. **Architecture mismatch**: Ensure 64-bit executable on 64-bit systems

## Testing

Test the portable executable by:
1. Copying the `exe\` folder to a different computer
2. Running `interfaz.exe` from that location
3. Verifying all functionality works correctly
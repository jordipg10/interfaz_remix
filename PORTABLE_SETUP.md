# Portable Executable Setup for `interfaz_remix.exe`

This document explains how to make the `interfaz_remix.exe` executable work on
other Windows computers that do not have a gfortran toolchain installed.

## The problem

`interfaz_remix.exe` is built with gfortran (MinGW-w64 / MSYS2 / TDM-GCC) and at
runtime depends on the following DLLs:

- `libgfortran_64-5.dll` — Fortran runtime
- `libgcc_s_seh_64-1.dll` — GCC runtime
- `libwinpthread_64-1.dll` — POSIX threads
- `libquadmath_64-0.dll`  — Quad-precision math

These DLLs are normally installed alongside gfortran but are not present on
machines without the toolchain.

## Repository layout

```
bin/   interfaz_remix.exe + the four runtime DLLs (this is what you ship)
lib/   Backup copy of the four runtime DLLs (used by copy_local_dlls.ps1)
```

The `bin/` folder in this repository already contains the executable together
with the four DLLs, so it is portable as-is.

## Solutions

### Option 1: Ship the executable together with the DLLs (recommended)

This is the simplest and recommended approach.

#### Automatic setup

- PowerShell (in-repo DLLs):
	```powershell
	./copy_local_dlls.ps1
	```
- PowerShell (locate DLLs from your gfortran installation):
	```powershell
	./copy_dlls.ps1
	```

Both scripts ensure the four DLLs are present next to `bin/interfaz_remix.exe`.

#### Manual setup

1. Locate the four DLLs. Preferred source: `lib/` in this repo. Otherwise look
   in your gfortran installation, e.g. `C:\TDM-GCC-64\bin\` or
   `C:\msys64\mingw64\bin\`.
2. Copy them next to the executable:
	```powershell
	Copy-Item lib\libgfortran_64-5.dll, lib\libgcc_s_seh_64-1.dll, `
			  lib\libwinpthread_64-1.dll, lib\libquadmath_64-0.dll `
			  bin\ -Force
	```

> The DLLs must live in the same directory as `interfaz_remix.exe` (or in a
> directory on `PATH`) for the program to start.

### Option 2: Static linking (optional)

You can try to embed the gfortran/GCC runtimes to avoid shipping DLLs. This is
not enabled by default. Suggested linker flags on Windows / MinGW:

```
-static-libgfortran -static-libgcc -static-libquadmath
```

Notes:

- Full `-static` is often problematic on Windows and is not recommended.
- Some library combinations may still pull in DLLs at runtime.
- To use this, modify the `link` task in `.vscode/tasks.json` (or add a
  separate `link-static` task).

### Option 3: Installer package

Wrap the executable plus its DLLs (and optionally `DB/`, `examples/`,
`documentation/`) in an installer:

- NSIS, InnoSetup, or WiX are all reasonable choices.
- The installer should place the four DLLs in the same folder as the executable.

## Distribution checklist

To distribute a portable build:

1. Build the executable (see `BUILD_GUIDE.md`, e.g. via the `rebuild` VS Code task).
2. Make sure `bin/` contains:
	- `interfaz_remix.exe`
	- `libgfortran_64-5.dll`
	- `libgcc_s_seh_64-1.dll`
	- `libwinpthread_64-1.dll`
	- `libquadmath_64-0.dll`
3. (Optional) Include the `DB/`, `examples/`, and `documentation/` folders so
   the user has databases and ready-to-run problems.
4. Zip the result and test it on a clean machine without gfortran installed.

## Helper scripts and tasks

PowerShell scripts at the repository root:

- `copy_local_dlls.ps1` — copies the four DLLs from `lib/` into `bin/`.
- `copy_dlls.ps1` — locates and copies the DLLs from your gfortran installation
  into `bin/`.
- `build_multiplatform.ps1` — multi-platform build driver (Docker-based for Linux).

VS Code tasks (see `.vscode/tasks.json`):

- `compile`, `compile-discr`, `compile-chem`, `compile-main` — build the objects.
- `link` — link `obj/*.o` into `bin/interfaz_remix.exe` (dynamic).
- `rebuild` — `clean` → `compile-discr` → `compile-chem` → `compile-main` → `link`.
- `run`, `link-and-run`, `compile-main-and-link-and-run` — launch the executable.
- `clean` — remove `obj/*.o` and `mod/*.mod`.

> The current `tasks.json` does not include dedicated `copy-dlls` /
> `build-with-dlls` tasks. Use the PowerShell scripts above instead.

## Troubleshooting

If the executable still does not run on another machine:

1. **Check Windows version / architecture**: ship the 64-bit build to 64-bit
   targets. Re-build for 32-bit if you need to support legacy systems.
2. **Verify dependencies**: from a Visual Studio Developer Prompt, run
   `dumpbin /dependents bin\interfaz_remix.exe` to list the DLLs the binary
   actually requires; make sure each one is present next to the executable
   (or on the system).
3. **Missing Visual C++ runtime**: some systems may also require the
   Microsoft Visual C++ Redistributable.
4. **Antivirus / SmartScreen**: an unsigned executable can be blocked on
   first run. Right-click → Properties → Unblock if needed.

## Testing portability

1. Copy the `bin/` folder (and any required data folders such as `DB/` and
   `examples/`) to a Windows machine without gfortran installed.
2. Open a terminal in that folder and run `.\interfaz_remix.exe`.
3. Walk through the interactive prompts described in `README.md` to confirm
   end-to-end functionality.

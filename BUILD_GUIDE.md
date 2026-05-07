# Cross-Platform Build Guide

This guide explains how to build the `interfaz_remix` executable on different
operating systems. The Fortran sources live under `src/`, compiled object files
go to `obj/`, generated module files to `mod/`, and the linked executable to
`bin/interfaz_remix.exe` (Windows) or `bin/interfaz_remix` (Linux / macOS).

## Prerequisites

- General
	- A recent gfortran toolchain installed and available on `PATH`.
	- VS Code (optional) if you want to use the provided tasks in `.vscode/tasks.json`.
	- Write access to create/use the folders: `obj/`, `mod/`, `bin/`.
- Windows
	- MinGW / MSYS2 / TDM-GCC gfortran (64-bit). The tasks in this repo currently
	  hard-code the path
	  `C:\Users\user2319\OneDrive\Documentos\fortran\mingw64\bin\gfortran.exe` — adjust
	  it to match your installation, or replace it with `gfortran` if it is on `PATH`.
	- Runtime DLLs needed next to the executable (already shipped in `bin/`):
		- `libgfortran_64-5.dll`
		- `libquadmath_64-0.dll`
		- `libgcc_s_seh_64-1.dll`
		- `libwinpthread_64-1.dll`
- Linux
	- Install `gfortran` via your distro's package manager.
- macOS
	- Install gfortran via Homebrew: `brew install gcc`.

## Repository layout (relevant for building)

```
src/   Fortran sources (.f90)
obj/   Compiled objects (.o) — build output
mod/   Generated Fortran modules (.mod) — build output
bin/   Linked executable + Windows runtime DLLs
lib/   Optional copy of the runtime DLLs used by the helper PowerShell scripts
```

## Windows

### Option 1: Use the VS Code tasks

The available tasks (see `.vscode/tasks.json`) are:

| Task                              | Purpose                                                                 |
|-----------------------------------|-------------------------------------------------------------------------|
| `create-obj-dir`                  | Create `obj/` if missing.                                               |
| `create-mod-dir`                  | Create `mod/` if missing.                                               |
| `compile-discr`                   | Compile the discretization / utilities layer.                           |
| `compile-chem`                    | Compile the chemistry layer.                                            |
| `compile-main`                    | Compile `src/main_interfaz.f90`.                                        |
| `compile`                         | Compile the full source tree in the correct order.                      |
| `compile-to-obj-dir`              | Compile a single file (template task) into `obj/`.                      |
| `compile-to-obj-dir-cwd`          | Compile a single file using `obj/` as cwd (template task).              |
| `link`                            | Link `obj/*.o` into `bin/interfaz_remix.exe` (dynamic).                 |
| `clean`                           | Remove `obj/*.o` and `mod/*.mod`.                                       |
| `rebuild`                         | `clean` → `compile-discr` → `compile-chem` → `compile-main` → `link`.   |
| `run`                             | Run `./bin/interfaz_remix.exe` (errors out if it does not exist).       |
| `debug`                           | Launch `gdb ./bin/interfaz_remix.exe`.                                  |
| `link-and-run`                    | `link` → `run`.                                                         |
| `compile-main-and-link-and-run`   | `compile-main` → `link` → `run`.                                        |

Typical flows:

- First build / full rebuild → run the `rebuild` task.
- Quick edit-compile-run cycle on the driver → run `compile-main-and-link-and-run`.
- After editing only the chemistry layer → run `compile-chem`, then `link-and-run`.

The Windows runtime DLLs are already present in `bin/`, so no extra copy step is
needed when running locally. If you move the executable elsewhere, copy the
DLLs along with it (see [Distribution](#distribution)).

### Option 2: Manual build (PowerShell)

```powershell
# 1. Prepare output directories
New-Item -ItemType Directory -Force -Path obj, mod, bin | Out-Null

# 2. Compile (module order matters; the VS Code `compile` task encodes the
#    full ordered list of source files. For a manual build, mirror that order
#    or invoke gfortran per file.)
cd obj
gfortran -g -c -O0 -fcheck=all -fbacktrace -J ../mod ../src/<file>.f90
# ... repeat for every source file in the correct dependency order ...

# 3. Link
gfortran -o ../bin/interfaz_remix.exe *.o

# 4. Ensure DLLs are next to the executable (already shipped in bin/, but if
#    you copied the exe somewhere else):
Copy-Item ..\lib\libgfortran_64-5.dll, ..\lib\libquadmath_64-0.dll, `
		  ..\lib\libgcc_s_seh_64-1.dll, ..\lib\libwinpthread_64-1.dll `
		  ..\bin\ -Force
```

The repository also ships two helper scripts:

- `copy_dlls.ps1` — locates the runtime DLLs from your gfortran installation
  and copies them next to the executable.
- `copy_local_dlls.ps1` — copies the DLLs from the in-repo `lib/` folder.

Use whichever matches your environment.

## Linux

```bash
# Install gfortran
sudo apt-get install gfortran          # Ubuntu/Debian
sudo yum install gcc-gfortran          # RHEL/CentOS
sudo dnf install gcc-gfortran          # Fedora

# Prepare output directories
mkdir -p obj mod bin

# Compile (mirror the order encoded in the VS Code `compile` task)
cd obj
gfortran -g -c -O0 -fcheck=all -fbacktrace -ffree-line-length-none \
		 -fno-range-check -J ../mod ../src/<file>.f90
# ... repeat for every source file in the correct dependency order ...

# Link (dynamic)
gfortran -o ../bin/interfaz_remix *.o

# Optional: static linking (may require extra system static libraries)
# gfortran -static -o ../bin/interfaz_remix *.o
```

A `Dockerfile.linux` is not currently included in the repo. If you add one,
the helper script `build_multiplatform.ps1` can be used to drive a Docker-based
multi-platform build.

## macOS

```bash
# Install gfortran via Homebrew
brew install gcc

# Prepare output directories
mkdir -p obj mod bin

# Compile (mirror the order encoded in the VS Code `compile` task)
cd obj
gfortran -g -c -O0 -fcheck=all -fbacktrace -ffree-line-length-none \
		 -fno-range-check -J ../mod ../src/<file>.f90
# ... repeat for every source file in the correct dependency order ...

# Link (dynamic)
gfortran -o ../bin/interfaz_remix *.o
```

Note: fully static linking is generally not supported on macOS.

## Distribution

### Windows

- Ship `bin/interfaz_remix.exe` together with the four runtime DLLs already in `bin/`:
	- `libgfortran_64-5.dll`
	- `libquadmath_64-0.dll`
	- `libgcc_s_seh_64-1.dll`
	- `libwinpthread_64-1.dll`
- Optionally include `DB/`, `documentation/`, and `examples/` so the user has
  databases and ready-to-run problems.
- Zip the `bin/` folder (plus the auxiliary folders above) for distribution.
- See `PORTABLE_SETUP.md` for additional portability notes.

### Linux / macOS

- Ship the single dynamically-linked executable. Make sure the standard
  gfortran / libc shared libraries exist on the target system.
- Set the executable bit: `chmod +x bin/interfaz_remix`.

## Available helper scripts

- `build_multiplatform.ps1` — multi-platform build driver (Docker-based for Linux).
- `copy_dlls.ps1` — copy gfortran runtime DLLs from your installation into `bin/`.
- `copy_local_dlls.ps1` — copy DLLs from the in-repo `lib/` folder into `bin/`.

## Troubleshooting

- **`gfortran` not found**
	- Ensure `gfortran` is on `PATH` (`gfortran --version`), or update the
	  hard-coded path in `.vscode/tasks.json` to match your installation.
- **VS Code task fails after refactoring sources**
	- Run `clean`, then `rebuild`. The order of files in the `compile` and
	  `link` tasks must match the module dependency graph.
- **Missing DLL at runtime (Windows)**
	- Verify the four DLLs listed above are next to `interfaz_remix.exe`.
	  If they are not, run `copy_local_dlls.ps1` (or `copy_dlls.ps1`).
- **Command-line length issues when invoking gfortran manually**
	- Compile in batches, or rely on the `compile` VS Code task which already
	  splits the work appropriately.
- **Spurious IEEE warnings on exit**
	- The driver clears IEEE flags before terminating; if you still see
	  warnings, check that your gfortran build supports `ieee_arithmetic` /
	  `ieee_exceptions`.

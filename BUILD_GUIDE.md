# Cross-Platform Build Guide

This guide explains how to build the `interfaz_remix` executable on different
operating systems.

| Folder | Role                                                      |
|--------|-----------------------------------------------------------|
| `src/` | Fortran sources (`.f90`).                                 |
| `obj/` | Compiled objects (`.o`) — build output.                   |
| `mod/` | Generated Fortran modules (`.mod`) — build output.        |
| `bin/` | Linked executable (+ Windows runtime DLLs).               |
| `lib/` | In-repo backup of the runtime DLLs (used by helper scripts). |

The linked executable is `bin/interfaz_remix.exe` on Windows or
`bin/interfaz_remix` on Linux / macOS.

## Prerequisites

### General

- A recent **gfortran** toolchain.
- **VS Code** (optional) — only needed to use the tasks in
  [.vscode/tasks.json](.vscode/tasks.json).
- Write access to create/use `obj/`, `mod/`, `bin/`.

### Windows

- MinGW-w64 / MSYS2 / TDM-GCC gfortran (64-bit). The VS Code tasks currently
  hard-code the path

  ```
  C:\Users\jordi\OneDrive\Documentos\fortran\mingw64\bin\gfortran.exe
  ```

  Edit the `command` field of every task in
  [.vscode/tasks.json](.vscode/tasks.json) if your gfortran lives elsewhere
  (for example `C:\TDM-GCC-64\bin\gfortran.exe` or
  `C:\msys64\mingw64\bin\gfortran.exe`), or replace it with plain `gfortran`
  if it is on `PATH`.

- Runtime DLLs needed next to the executable (already shipped in `bin/`):
  - `libgfortran_64-5.dll`
  - `libquadmath_64-0.dll`
  - `libgcc_s_seh_64-1.dll`
  - `libwinpthread_64-1.dll`

### Linux

Install gfortran via your distro's package manager:

```bash
sudo apt-get install gfortran          # Ubuntu/Debian
sudo yum install gcc-gfortran          # RHEL/CentOS
sudo dnf install gcc-gfortran          # Fedora
```

### macOS

```bash
brew install gcc
```

## Windows

### Option 1 — VS Code tasks

The available tasks (defined in [.vscode/tasks.json](.vscode/tasks.json)):

| Task                       | Purpose                                                                                  |
|----------------------------|------------------------------------------------------------------------------------------|
| `create-obj-dir`           | Create `obj/` if missing.                                                                |
| `create-mod-dir`           | Create `mod/` if missing.                                                                |
| `compile-discr`            | Compile the discretization / utilities layer.                                            |
| `compile-chem`             | Compile the chemistry layer.                                                             |
| `compile-main`             | Compile `src/main_interfaz.f90`.                                                         |
| `compile`                  | Compile the full source tree in the correct order.                                       |
| `compile-to-obj-dir`       | Template task: compile a single file into `obj/`.                                        |
| `compile-to-obj-dir-cwd`   | Template task: compile a single file using `obj/` as cwd.                                |
| `link`                     | Link `obj/*.o` into `bin/interfaz_remix.exe` (dynamic).                                  |
| `clean`                    | Remove `obj/*.o` and `mod/*.mod`.                                                        |
| `rebuild`                  | `clean` → `compile-discr` → `compile-chem` → `compile-main` → `link`.                    |
| `run`                      | Run `./bin/interfaz_remix.exe` (errors out if the executable is missing).                |
| `debug`                    | Launch `gdb ./bin/interfaz_remix.exe`.                                                   |
| `link-and-run`             | `link` → `run`.                                                                          |
| `compile-main-and-link`    | `compile-main` → `link`.                                                                 |

Typical flows:

- First build / full rebuild → run the **`rebuild`** task.
- After editing only the driver → run **`compile-main-and-link`**, then `run`.
- After editing only the chemistry layer → run **`compile-chem`**, then
  `link-and-run`.

The Windows runtime DLLs are already present in `bin/`, so no extra copy step
is needed when running locally. If you move the executable elsewhere, copy the
DLLs along with it (see [Distribution](#distribution)).

### Option 2 — Manual build (PowerShell)

```powershell
# 1. Prepare output directories
New-Item -ItemType Directory -Force -Path obj, mod, bin | Out-Null

# 2. Compile (module order matters — mirror the order encoded in the
#    `compile-discr`, `compile-chem`, `compile-main` tasks of
#    .vscode/tasks.json, or call the `compile` task which does it all).
cd obj
gfortran -g -c -O0 -fcheck=all -fbacktrace -J ../mod ../src/<file>.f90
# ... repeat for every source file in the correct dependency order ...

# 3. Link
gfortran -o ../bin/interfaz_remix.exe *.o

# 4. Make sure the four DLLs are next to the executable (already shipped in
#    bin/, but if you copied the exe somewhere else):
Copy-Item ..\lib\libgfortran_64-5.dll, ..\lib\libquadmath_64-0.dll, `
          ..\lib\libgcc_s_seh_64-1.dll, ..\lib\libwinpthread_64-1.dll `
          ..\bin\ -Force
```

Helper scripts:

- `copy_dlls.ps1` — locates the runtime DLLs in your gfortran installation and
  copies them next to the executable.
- `copy_local_dlls.ps1` — copies the DLLs from the in-repo `lib/` folder.

## Linux

```bash
mkdir -p obj mod bin

cd obj
gfortran -g -c -O0 -fcheck=all -fbacktrace -ffree-line-length-none \
         -fno-range-check -J ../mod ../src/<file>.f90
# ... repeat for every source file in the correct dependency order ...

# Link (dynamic)
gfortran -o ../bin/interfaz_remix *.o

# Optional: static linking (may require extra system static libraries)
# gfortran -static -o ../bin/interfaz_remix *.o
```

`build_multiplatform.ps1` is provided as a Docker-based driver for
multi-platform builds. A `Dockerfile.linux` is not currently included; add
one if you want to drive the Linux build through Docker.

## macOS

```bash
brew install gcc

mkdir -p obj mod bin

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

- Ship `bin/interfaz_remix.exe` together with the four runtime DLLs already in
  `bin/`:
  - `libgfortran_64-5.dll`
  - `libquadmath_64-0.dll`
  - `libgcc_s_seh_64-1.dll`
  - `libwinpthread_64-1.dll`
- Optionally include `DB/`, `documentation/`, and `examples/` so the user has
  databases and ready-to-run problems.
- Zip the `bin/` folder (plus the auxiliary folders above) for distribution.
- See [PORTABLE_SETUP.md](PORTABLE_SETUP.md) for additional portability notes.

### Linux / macOS

- Ship the dynamically-linked executable. Make sure the standard
  gfortran / libc shared libraries exist on the target system.
- Set the executable bit: `chmod +x bin/interfaz_remix`.

## Helper scripts

- `build_multiplatform.ps1` — multi-platform build driver (Docker-based for Linux).
- `copy_dlls.ps1` — copy gfortran runtime DLLs from your installation into `bin/`.
- `copy_local_dlls.ps1` — copy DLLs from the in-repo `lib/` folder into `bin/`.

## Troubleshooting

- **`gfortran` not found**
  Ensure `gfortran` is on `PATH` (`gfortran --version`), or update the
  hard-coded path in [.vscode/tasks.json](.vscode/tasks.json) to match your
  installation.
- **`Cannot open module file 'xxx.mod'` during compilation**
  The compile order is wrong or a source file is missing from the relevant
  `compile-*` task. Run `clean`, then `rebuild`. The order of files in the
  `compile-*` and `link` tasks must match the module dependency graph.
- **VS Code task fails after refactoring sources**
  Same as above — run `clean`, then `rebuild`. If you renamed a source file,
  remember to update the corresponding entries in `compile-*` (source name)
  and `link` (object name) inside [.vscode/tasks.json](.vscode/tasks.json).
- **Missing DLL at runtime (Windows)**
  Verify the four DLLs listed above are next to `interfaz_remix.exe`. If they
  are not, run `copy_local_dlls.ps1` (or `copy_dlls.ps1`).
- **Command-line length issues when invoking gfortran manually**
  Compile in batches, or rely on the VS Code `compile` task which already
  splits the work appropriately.
- **Spurious IEEE warnings on exit**
  The driver clears IEEE flags before terminating; if you still see warnings,
  check that your gfortran build supports `ieee_arithmetic` /
  `ieee_exceptions`.

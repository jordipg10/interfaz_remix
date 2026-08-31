# Cross-Platform Build Guide

This guide explains how to build the `interfaz_remix` executable on different
operating systems.

| Folder | Role                                                      |
|--------|-----------------------------------------------------------|
| `src/` | Fortran sources (`.f90`).                                 |
| `obj/` | Compiled objects (`.o`) — build output.                   |
| `mod/` | Generated Fortran modules (`.mod`) — build output.        |
| `bin/` | Linked executable (standalone — no extra DLLs required).  |
| `lib/` | Static LAPACK/BLAS archives used at link time.            |

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

- No runtime DLLs are required next to the executable. The `link` VS Code
  task passes `-static`, which embeds the gfortran, GCC, quadmath, and
  winpthread runtimes directly into `bin/interfaz_remix.exe`. The only
  remaining dependencies are standard Windows system DLLs (`KERNEL32.dll`
  and the Universal CRT) that are always present on Windows 10/11.

### Linux

Install gfortran and the LAPACK/BLAS development packages via your distro's
package manager (the bundled `lib/*.a` archives are Windows-only and are NOT
used on Linux):

```bash
sudo apt-get install gfortran liblapack-dev libblas-dev   # Ubuntu/Debian
sudo yum install gcc-gfortran lapack-devel blas-devel       # RHEL/CentOS
sudo dnf install gcc-gfortran lapack-devel blas-devel       # Fedora
```

### macOS

```bash
brew install gcc   # provides gfortran; LAPACK/BLAS come from the built-in Accelerate framework
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
| `link`                     | Link `obj/*.o` into `bin/interfaz_remix.exe` (static — self-contained).                  |
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

The executable is statically linked, so no extra DLLs are needed when running
locally or when copying `bin/interfaz_remix.exe` to another machine.

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

# 3. Link (static — no runtime DLLs needed on the target machine)
gfortran -static -o ../bin/interfaz_remix.exe *.o `
         ${workspaceFolder}/lib/liblapack.a `
         ${workspaceFolder}/lib/librefblas.a
```

Helper scripts:

- `copy_dlls.ps1` / `copy_local_dlls.ps1` — legacy scripts for copying
  runtime DLLs. Not needed when building with `-static` (the default).

## Linux

### Option 1 - build.sh (recommended)

From the repository root:

```bash
sudo apt-get install -y gfortran liblapack-dev libblas-dev   # once
chmod +x build.sh
./build.sh                                # -> bin/interfaz_remix
FFLAGS="-O2 -fbacktrace" ./build.sh       # optimised/release build
```

`build.sh` compiles every source in the correct module order (mirroring the
`compile-*` VS Code tasks) and links against the system LAPACK/BLAS. From a
Windows PC, run it inside WSL (Ubuntu).

### Option 2 - manual

```bash
mkdir -p obj mod bin

cd obj
gfortran -g -c -O0 -fcheck=all -fbacktrace -ffree-line-length-none \
         -fno-range-check -J ../mod ../src/<file>.f90
# ... repeat for every source file in the correct dependency order ...

# Link against the SYSTEM LAPACK/BLAS (NOT the Windows lib/*.a archives).
# -static-libgfortran embeds the gfortran runtime; LAPACK/BLAS stay dynamic.
gfortran -static-libgfortran -static-libgcc -o ../bin/interfaz_remix *.o -llapack -lblas
```

## macOS

### Option 1 - build.sh (recommended)

```bash
brew install gcc          # provides gfortran; Accelerate is built in
chmod +x build.sh
./build.sh                # -> bin/interfaz_remix
```

### Option 2 - manual

```bash
mkdir -p obj mod bin

cd obj
gfortran -g -c -O0 -fcheck=all -fbacktrace -ffree-line-length-none \
         -fno-range-check -J ../mod ../src/<file>.f90
# ... repeat for every source file in the correct dependency order ...

# Link against LAPACK/BLAS from the built-in Accelerate framework.
gfortran -o ../bin/interfaz_remix *.o -framework Accelerate
```

Note: fully static linking is not supported on macOS. Distribute the dynamic
binary alongside a note that gfortran (Homebrew `gcc`) must be installed.

## Continuous Integration (Linux + macOS binaries)

The repository includes a GitHub Actions workflow,
[.github/workflows/build.yml](.github/workflows/build.yml), that builds the
executable natively on `ubuntu-latest` and `macos-latest` and uploads each as a
downloadable artifact. This is the simplest way to obtain Linux and macOS
binaries from a Windows-only machine:

1. Push the repository to GitHub, or open the **Actions** tab and choose
   **build -> Run workflow**.
2. When the run finishes, download `interfaz_remix-ubuntu-latest` and
   `interfaz_remix-macos-latest` from the run's **Artifacts** section.

## Distribution

### Windows

- Ship **only** `bin/interfaz_remix.exe` — no DLLs required. The executable
  is statically linked and carries all gfortran runtimes internally.
- Optionally include `DB/`, `documentation/`, and `examples/` so the user has
  databases and ready-to-run problems.
- Zip the `bin/` folder (plus the auxiliary folders above) for distribution.
- See [PORTABLE_SETUP.md](PORTABLE_SETUP.md) for additional portability notes.
- The only system requirement on the target PC is 64-bit Windows 10/11 (or
  Windows 7/8.1 with the Universal C Runtime update applied).

### Linux / macOS

- Ship the dynamically-linked executable. On Linux the target system needs the
  LAPACK/BLAS shared libraries (`liblapack.so.3`, `libblas.so.3`) plus the
  standard libc; on macOS LAPACK/BLAS come from the always-present Accelerate
  framework. (`build.sh` embeds the gfortran runtime via `-static-libgfortran`
  on Linux.)
- Set the executable bit: `chmod +x bin/interfaz_remix`.

## Running on another PC

The executable is **fully self-contained**. The target PC does **not** need
gfortran (or any Fortran compiler) installed — all gfortran runtimes are
statically embedded into `bin/interfaz_remix.exe` by the `-static` link flag.

What the target PC **does** need:

- 64-bit Windows 10 or 11 (or Windows 7/8.1 with the Universal CRT update
  — KB2999226 — applied).
- Only `interfaz_remix.exe` itself; no companion DLLs are required.

The `DB/` directory (or the path configured inside the program) must be
accessible at runtime so the program can open its database files.

## Helper scripts

- `build_multiplatform.ps1` — packages `bin/interfaz_remix.exe` into `dist/`
  for Windows, and (if Docker plus a `Dockerfile.linux` are available) builds
  a Linux executable in the same `dist/` folder. It does **not** rebuild the
  Windows binary itself — run the `rebuild` VS Code task first.
- `copy_dlls.ps1` / `copy_local_dlls.ps1` — legacy scripts for copying
  gfortran runtime DLLs. Not needed when linking with `-static` (the default).

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
  This should not happen with the statically-linked build. If you see this
  error, the exe was likely built without `-static` by a different toolchain.
  Rebuild using the `rebuild` VS Code task, which passes `-static`.
- **Command-line length issues when invoking gfortran manually**
  Compile in batches, or rely on the VS Code `compile` task which already
  splits the work appropriately.
- **Spurious IEEE warnings on exit**
  The driver clears IEEE flags before terminating; if you still see warnings,
  check that your gfortran build supports `ieee_arithmetic` /
  `ieee_exceptions`.

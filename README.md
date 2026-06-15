# interfaz_remix

`interfaz_remix` provides the WMA (Water Mixing Approach) interfaces that can
be called from a reactive transport simulation. The intended workflow is:

1. The user solves the conservative transport step externally.
2. The resulting concentrations are written to a file.
3. `interfaz_remix` is invoked at every time step to apply the reactive mixing
   stage on top of those concentrations.

## Repository layout

| Folder / file               | Purpose                                                                 |
|-----------------------------|-------------------------------------------------------------------------|
| `src/`                      | Fortran source code (`.f90`).                                           |
| `obj/`                      | Compiled object files (build output).                                   |
| `mod/`                      | Generated Fortran module files (build output).                          |
| `bin/`                      | Linked Windows executable (standalone — no extra DLLs required).        |
| `lib/`                      | Static LAPACK/BLAS archives used at link time.                           |
| `DB/`                       | Chemical databases used by the program.                                 |
| `examples/`                 | Sample reactive-transport problems.                                     |
| `documentation/`            | Documentation for the main program, the interfaces and the input files. |
| `.vscode/`                  | VS Code build/run/debug tasks.                                          |
| `BUILD_GUIDE.md`            | Cross-platform build instructions.                                      |
| `PORTABLE_SETUP.md`         | How to ship the executable to machines without gfortran.                |
| `build_multiplatform.ps1`   | Helper script for multi-platform builds.                                |
| `copy_dlls.ps1`             | Legacy: copy gfortran runtime DLLs from the toolchain (not needed with static build). |
| `copy_local_dlls.ps1`       | Legacy: copy in-repo runtime DLLs (not needed with static build).       |

## Quick start (Windows)

The prebuilt executable lives at `bin/interfaz_remix.exe`. It is statically
linked — no companion DLLs or gfortran installation are required on the
running machine.

To run it:

```powershell
cd .\bin
.\interfaz_remix.exe
```

You can also launch it from VS Code with the `run` task defined in
[.vscode/tasks.json](.vscode/tasks.json).

## Prerequisites (for building from source)

- A recent **gfortran** toolchain (64-bit on Windows).
- Optional: **VS Code**, to use the tasks shipped in `.vscode/tasks.json`.
- Write access to create/use the folders `obj/`, `mod/`, and `bin/`.

### gfortran path on Windows

The VS Code tasks currently hard-code the path

```
C:\Users\jordi\OneDrive\Documentos\fortran\mingw64\bin\gfortran.exe
```

If your gfortran lives somewhere else (for example `C:\TDM-GCC-64\bin\` or
`C:\msys64\mingw64\bin\`), edit the `command` field of every task in
[.vscode/tasks.json](.vscode/tasks.json) accordingly, or replace it with plain
`gfortran` if the compiler is on your `PATH`.

On Linux/macOS, install gfortran through your package manager
(`apt`, `dnf`, `brew`, …) and make sure it is on `PATH`.

## Building

Detailed cross-platform build instructions live in
[BUILD_GUIDE.md](BUILD_GUIDE.md). The most useful VS Code tasks are:

| Task                     | Purpose                                                                 |
|--------------------------|-------------------------------------------------------------------------|
| `create-obj-dir`         | Create the `obj/` folder if it does not exist.                          |
| `create-mod-dir`         | Create the `mod/` folder if it does not exist.                          |
| `compile-discr`          | Compile the discretization / utilities layer.                           |
| `compile-chem`           | Compile the chemistry layer.                                            |
| `compile-main`           | Compile `src/main_interfaz.f90`.                                        |
| `compile`                | Compile the full source tree in the correct order.                      |
| `link`                   | Link `obj/*.o` into `bin/interfaz_remix.exe`.                           |
| `clean`                  | Remove `obj/*.o` and `mod/*.mod`.                                       |
| `rebuild`                | `clean` → `compile-discr` → `compile-chem` → `compile-main` → `link`.   |
| `run`                    | Run `./bin/interfaz_remix.exe`.                                         |
| `debug`                  | Launch `gdb` on `./bin/interfaz_remix.exe`.                             |
| `link-and-run`           | `link` → `run`.                                                         |
| `compile-main-and-link`  | `compile-main` → `link`.                                                |

## Distributing the executable

Because the `link` task uses `-static`, all gfortran runtimes are baked into
`bin/interfaz_remix.exe`. To distribute, copy just the executable (and the
`DB/` database folder). No companion DLLs are needed. See
[PORTABLE_SETUP.md](PORTABLE_SETUP.md) and [BUILD_GUIDE.md](BUILD_GUIDE.md)
for full details.

## Examples

The `examples/` folder contains small reactive-transport test cases that can be
freely modified:

- `calcite_eq/`
- `cc_anh_eq/`
- `denit_2reacts/`
- `denit_2reacts_calcite/`
- `denit_ext/`
- `gypsum_eq/`
- `gypsum_kin/`

## Documentation

The `documentation/` folder describes the main program (`main_interfaz.f90`),
the public interfaces (`interfaz_comps_arch`, `interfaz_esp_arch`) and the
required input files. The source code is annotated with Doxygen-style comments
(`!>`, `!<`), so HTML/LaTeX documentation can be generated with
[Doxygen](https://www.doxygen.nl/).

The `DB/` folder contains the chemical databases used by the program.

## Execution flow (summary)

The driver in [src/main_interfaz.f90](src/main_interfaz.f90) runs **one**
reactive-mixing iteration. The prompts (in order) are:

1. Launch the executable (`./bin/interfaz_remix.exe`).
2. Database path: directory containing the chemical database.
3. Problem path: directory containing the problem-specific input files.
4. Root of the input and output files (filename prefix used by the driver).
5. Number of targets in the mesh.
6. Whether the file with the aqueous component concentrations of the initial
   and external water types has already been generated (`1`: yes, `0`: no).
   - If `0`, provide the name of the file to **write** those concentrations to.
   - If `1`, provide the name of the **existing** file to read them from.
7. Name of the file where post-mixing concentrations should be written.
8. Name of the file (inside the problem directory) that contains `u_tilde` —
   the concentrations after one conservative-transport step.
9. Layout of the `u_tilde` file: `1` if rows are targets and columns are
   components, `0` if rows are components and columns are targets.
10. Time step (`Δt > 0`).

The program then runs one reactive-mixing iteration and exits. To drive
multiple time steps, invoke `interfaz_remix.exe` once per step from the
outer transport loop (refreshing the `u_tilde` file between calls).

Input can also be redirected from a file (e.g. `interfaz_remix.exe < fort.5`).
Blank lines and full-line comments starting with `!` or `#` are ignored, so
the input file can be annotated.

The program automatically picks the right reactive-mixing interface based on
the chemical system, and on the orientation flag of step 9 when applicable:

- `interfaz_esp_arch` — no equilibrium reactions (kinetic only).
- `interfaz_comps_arch_eq` / `interfaz_comps_arch_eq_T` — equilibrium
  reactions only (no kinetics). The `_T` variant is used when the input
  matrix is transposed (rows = targets).
- `interfaz_comps_arch_eq_kin` / `interfaz_comps_arch_eq_kin_T` — both
  equilibrium and kinetic reactions. The `_T` variant is used when the input
  matrix is transposed (rows = targets).

## Notes

- **Important:** when typing directories at the prompt, terminate them with a
  backslash (`\`) on Windows or a forward slash (`/`) on Linux/macOS.
- The program configures IEEE exceptions on start-up and clears them on exit
  to avoid spurious messages from the Fortran runtime.
- In non-interactive environments (CI, redirected stdin) the `read` calls
  report EOF and the program terminates with a clear message instead of
  crashing.

## Troubleshooting

- **`gfortran` not found / wrong path in VS Code tasks**
  Update the `command` field of every task in
  [.vscode/tasks.json](.vscode/tasks.json) to point at your gfortran, or use
  plain `gfortran` if it is on `PATH`.
- **`Cannot open module file 'xxx.mod'` during compilation**
  The compile order is wrong or a source file is missing from the relevant
  `compile-*` task. Run `clean`, then `rebuild`. The order of files in the
  `compile-*` and `link` tasks must match the module dependency graph.
- **Missing DLL at runtime (Windows)**
  This should not occur with the statically-linked build. If it does, the
  executable was likely built without `-static`. Rebuild using the `rebuild`
  VS Code task, which passes `-static` at link time.
- **Spurious IEEE warnings on exit**
  The driver clears IEEE flags before terminating; if you still see warnings,
  check that your gfortran build supports `ieee_arithmetic` /
  `ieee_exceptions`.

## License

This repository does not yet ship a `LICENSE` file. Until one is added, the
source code is to be treated as **"all rights reserved"**: please contact the
author before redistributing, reusing or publishing derivative work based on
it.

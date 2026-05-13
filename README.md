# interfaz_remix

`interfaz_remix` provides the WMA (Water Mixing Approach) interfaces that can
be called from a reactive transport simulation. The intended workflow is:

1. The user solves the conservative transport step externally.
2. The resulting concentrations are written to a file.
3. `interfaz_remix` is invoked at every time step to apply the reactive mixing
   stage on top of those concentrations.

## Repository layout

| Folder           | Purpose                                                                 |
|------------------|-------------------------------------------------------------------------|
| `src/`           | Fortran source code (`.f90`).                                           |
| `obj/`           | Compiled object files (build output).                                   |
| `mod/`           | Generated Fortran module files (build output).                          |
| `bin/`           | Linked Windows executable plus the runtime DLLs.                        |
| `lib/`           | Backup copy of the Windows runtime DLLs.                                |
| `DB/`            | Chemical databases used by the program.                                 |
| `examples/`      | Sample reactive-transport problems.                                     |
| `documentation/` | Documentation for the main program, the interfaces and the input files. |
| `.vscode/`       | VS Code build/run tasks.                                                |

## Quick start (Windows)

The prebuilt executable lives at `bin/interfaz_remix.exe` and is shipped next
to the four runtime DLLs it needs (`libgcc_s_seh_64-1.dll`,
`libgfortran_64-5.dll`, `libquadmath_64-0.dll`, `libwinpthread_64-1.dll`).

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

| Task                              | Purpose                                                                 |
|-----------------------------------|-------------------------------------------------------------------------|
| `compile-discr`                   | Compile the discretization / utilities layer.                           |
| `compile-chem`                    | Compile the chemistry layer.                                            |
| `compile-main`                    | Compile `src/main_interfaz.f90`.                                        |
| `compile`                         | Compile the full source tree in the correct order.                      |
| `link`                            | Link `obj/*.o` into `bin/interfaz_remix.exe`.                           |
| `clean`                           | Remove `obj/*.o` and `mod/*.mod`.                                       |
| `rebuild`                         | `clean` → `compile-discr` → `compile-chem` → `compile-main` → `link`.   |
| `run`                             | Run `./bin/interfaz_remix.exe`.                                         |
| `link-and-run`                    | `link` → `run`.                                                         |
| `compile-main-and-link-and-run`   | `compile-main` → `link` → `run`.                                        |

## Distributing the executable

Instructions for shipping the executable to machines that do not have gfortran
installed are in [PORTABLE_SETUP.md](PORTABLE_SETUP.md). The helper scripts
`copy_dlls.ps1` and `copy_local_dlls.ps1` automate copying the four runtime
DLLs next to the executable, and `build_multiplatform.ps1` covers
multi-platform builds.

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

1. Launch the executable (`./bin/interfaz_remix.exe`).
2. Provide the database directory, the problem directory and the `root` for
   the input/output files.
3. Enter the number of targets in the mesh.
4. Provide the file where the initial and external water-type component
   concentrations should be written.
5. Provide the file where the post-mixing concentrations should be written.
6. Provide the file (inside the problem directory) that contains `u_tilde` —
   the concentrations after one conservative-transport step.
7. Specify the layout of that file: enter `1` if rows are targets and columns
   are components, or `0` if rows are components and columns are targets.
8. Enter the initial time step (`Δt > 0`) and choose whether it is constant
   (`1`) or variable (`0`).
9. At every iteration, refresh the `u_tilde` file with the new transport
   solution and answer `1` to continue or `0` to stop. With a variable time
   step, the new `Δt` is requested at every iteration.

The program automatically picks the right reactive-mixing interface based on
the chemical system, and on the orientation flag of step 7 when applicable:

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
  Make sure the four runtime DLLs are next to `interfaz_remix.exe`. If they
  are not, run `copy_local_dlls.ps1` (in-repo DLLs) or `copy_dlls.ps1`
  (DLLs from your gfortran installation).
- **Spurious IEEE warnings on exit**
  The driver clears IEEE flags before terminating; if you still see warnings,
  check that your gfortran build supports `ieee_arithmetic` /
  `ieee_exceptions`.

## License

No license file is currently included in this repository. Until one is added,
treat the source as "all rights reserved" and contact the author before
redistributing or reusing it.

# Portable Executable Setup for `interfaz_remix.exe`

This document explains how `interfaz_remix.exe` is made fully self-contained
so it runs on any Windows 10/11 machine without a gfortran toolchain installed.

## How it works

The `link` VS Code task passes `-static` to gfortran at link time. This
statically embeds the following runtimes **inside** the executable itself:

- gfortran runtime (`libgfortran`)
- GCC runtime (`libgcc`)
- POSIX threads (`libwinpthread`)
- Quad-precision math (`libquadmath`)

The resulting `bin/interfaz_remix.exe` therefore has **no external DLL
dependencies** beyond the standard Windows system libraries (`KERNEL32.dll`
and the Universal CRT), which are always present on Windows 10/11.

You can verify this with:

```powershell
objdump -p bin\interfaz_remix.exe | Select-String "DLL Name:"
```

Expected output (KERNEL32 + Universal CRT only — no libgfortran / libgcc):

```
DLL Name: KERNEL32.dll
DLL Name: api-ms-win-crt-*.dll
...
```

## Distributing the executable

1. Build the executable (see [BUILD_GUIDE.md](BUILD_GUIDE.md), e.g. via the
   `rebuild` VS Code task).
2. Copy **only** `bin/interfaz_remix.exe` to the target machine — no companion
   DLLs are required.
3. (Optional) Include the `DB/`, `examples/`, and `documentation/` folders so
   the user has databases and ready-to-run problems.
4. Test on a clean machine without gfortran installed.

## System requirements on the target PC

- 64-bit Windows 10 or 11.
- Windows 7/8.1 also works provided the Universal CRT update (KB2999226) has
  been applied.
- No Fortran compiler or runtime library installation is needed.

## Installer package (optional)

If you prefer a polished distribution, wrap `interfaz_remix.exe` (plus `DB/`,
`examples/`, `documentation/`) in an installer built with NSIS, InnoSetup, or
WiX. No DLL packaging step is required.

## Rebuilding with static linking

The `-static` flag is already present in the `link` task in
`.vscode/tasks.json`. To rebuild:

```powershell
# From the repository root (VS Code task)
# Run the "rebuild" task, or manually:

cd obj
gfortran -static -o ..\bin\interfaz_remix.exe *.o `
         ..\lib\liblapack.a ..\lib\librefblas.a
```

On **Linux**, static linking against glibc is not recommended (and rarely
needed). Build natively; the resulting binary will depend on the system glibc
which is always present. On **macOS**, fully static linking is not supported —
distribute the dynamic binary.

## Legacy helper scripts

`copy_dlls.ps1` and `copy_local_dlls.ps1` were used when the executable was
built dynamically. They are no longer needed with the current static build but
are kept for reference.

## Troubleshooting

- **DLL error on startup** — The executable was linked without `-static` (by
  a different build). Rebuild using the `rebuild` VS Code task.
- **"The application was unable to start correctly (0xc0150002)"** — Universal
  CRT is missing. Apply Windows Update KB2999226.
- **Antivirus / SmartScreen** — An unsigned executable may be blocked on first
  run. Right-click → Properties → Unblock if needed.
- **Verify dependencies** — Use `objdump -p bin\interfaz_remix.exe | Select-String "DLL Name:"` to list the actual DLL imports.

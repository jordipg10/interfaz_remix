En este repositorio se encuentran las interfaces del WMA (Water Mixing Approach) para ser llamadas durante una simulación de transporte reactivo. La idea es que el usuario resuelva el transporte conservativo, escriba las concentraciones obtenidas en un archivo y después ejecute este programa en cada iteración temporal para aplicar la mezcla reactiva.

El ejecutable para Windows se encuentra en la carpeta `bin` y se llama `interfaz_remix.exe`. Está junto con las DLLs de gfortran (`libgcc_s_seh_64-1.dll`, `libgfortran_64-5.dll`, `libquadmath_64-0.dll`, `libwinpthread_64-1.dll`) que necesita durante la ejecución. En Windows, para ejecutarlo hay que abrir una terminal, ir al directorio `.\bin` y escribir `.\interfaz_remix.exe`. También se puede lanzar desde VS Code mediante la tarea `run` definida en `.vscode/tasks.json`.

En el caso que el usuario quiera compilarlo y enlazarlo, puede seguir las instrucciones que se encuentran en el archivo `BUILD_GUIDE.md`. Ahí se explican diferentes formas de compilar el programa en distintos sistemas operativos (Windows, Linux, macOS). Las tareas de VS Code (`compile-discr`, `compile-chem`, `compile-main`, `link`, `rebuild`, `run`, `clean`, `link-and-run`, `compile-main-and-link-and-run`) automatizan los pasos más habituales.

Las instrucciones para hacer el ejecutable transportable a otros ordenadores están en el archivo `PORTABLE_SETUP.md`. Si el usuario no tiene instalado gfortran, tendrá que copiar las DLLs necesarias junto al ejecutable para que se enlacen dinámicamente. Los scripts de ayuda `copy_dlls.ps1` y `copy_local_dlls.ps1` automatizan esa copia, y `build_multiplatform.ps1` cubre la compilación multiplataforma.

En la carpeta `examples` se encuentran varios ejemplos sencillos de transporte reactivo (`calcite_eq`, `cc_anh_eq`, `denit_2reacts`, `denit_ext`, `gypsum_eq`, `gypsum_kin`). Los archivos se pueden modificar tanto como el usuario desee.

En la carpeta `documentation` se encuentra la documentación del programa principal (`main_interfaz.f90`), de las interfaces (`interfaz_comps_arch`, `interfaz_esp_arch`) y de los archivos de entrada necesarios. El código fuente está anotado con comentarios estilo Doxygen (`!>`, `!<`), por lo que es posible generar documentación HTML/LaTeX con [Doxygen](https://www.doxygen.nl/).

El directorio `DB` contiene las bases de datos químicas que se usan en este programa.

## Flujo de ejecución resumido

1. Lanzar el ejecutable (`./bin/interfaz_remix.exe`).
2. Indicar el directorio de la base de datos, el directorio del problema y el `root` de los ficheros de entrada/salida.
3. Indicar el fichero donde escribir las concentraciones de los tipos de agua iniciales y externas.
4. Indicar el fichero donde escribir las concentraciones después de la mezcla reactiva.
5. Indicar el fichero (en el directorio del problema) que contiene `u_tilde` —las concentraciones después de un paso de transporte conservativo. Debe tener tantas filas como componentes y tantas columnas como targets.
6. Introducir el paso de tiempo inicial (`Δt > 0`) y elegir si será constante (`1`) o variable (`0`).
7. En cada iteración, actualizar el fichero de `u_tilde` con la nueva solución de transporte y responder `1` para continuar o `0` para terminar. Si el paso es variable, se pedirá el nuevo `Δt` en cada iteración.

El programa selecciona automáticamente la interfaz adecuada según el sistema químico:

- `interfaz_esp_arch` cuando no hay reacciones en equilibrio.
- `interfaz_comps_arch` cuando hay reacciones en equilibrio.

## Notas

- IMPORTANTE: a la hora de escribir los directorios en el terminal, hay que incluir un backslash (`\`) al final para usuarios de Windows o un frontslash (`/`) para usuarios de Linux y macOS.
- El programa configura las excepciones IEEE al arrancar y las limpia al terminar para evitar mensajes espurios del runtime de Fortran.
- En entornos no interactivos (CI, redirección de stdin) los `read` reportan EOF y el programa termina con un mensaje claro en vez de hacer crash.



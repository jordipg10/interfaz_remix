En este repositorio se encuentran las interfaces del WMA para ser llamadas durante una simulación de transporte reactivo. La idea es que el usuario resuelva el transporte conservativo, escriba las concentraciones obtenidas en un archivo y después ejecute este programa en cada iteración temporal.

El ejecutable para Windows se encuentra en la carpeta "exe" y se llama "interfaz.exe". Está junto con las librerías de gfortran que necesita durante la ejecución. En Windows, para ejecutarlo hay que entrar en el terminal, ir al directorio "..\exe" y escribir ".\interfaz.exe".

En el caso que el usuario quiera compilarlo y enlazarlo, puede seguir las instrucciones que se encuentran en el archivo "BUILD_GUIDE.md". Ahí se explican diferentes formas de compilar el programa en diferentes sistemas operativos (Windows, Linux, MacOS).

Las instrucciones de ejecución están en el archivo "PORTABLE_SETUP.exe" en el que se explica cómo hacer este programa transportable a otros ordenadores. Si el usuario no tiene instalado el gfortran, tendrá que copiar las DLLs necesarias en la carpeta "exe" para que sean enlazadas dinámicamente.

En la carpeta "examples" se encuentran varios ejemplos sencillos de transporte reactivo. Los archivos se pueden modificar tanto como el usuario desee.

En la carpeta "documentation" se encuentra la documentación del "main", de las interfaces y de los archivos de entrada necesarios.

El directorio "DB" contiene las bases de datos químicas que se usan en este programa.

IMPORTANTE: a la hora de escribir los directorios en el terminal, hay que incluir un backslash al final para usuarios de Windows o un frontslash para usuarios de Linux y MacOS.


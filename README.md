# SARS-CoV-2 y Minería de Datos

Taller para el colectivo [R Ladies Guadalajara](https://www.meetup.com/es-ES/rladies-guadalajara/events/279129632/) donde se explica un poco sobre cómo realizo Minería de Datos en nuestra investigación sobre la proteína Spike del virus SARS-CoV-2. En el repositorio se podrán encontrar los siguientes archivos:
- `Secuencias.zip`: Archivo comprimido con las secuencias en formato fasta
- `Secuencias.csv`: Archivo con los metadatos de las secuencias
- `SARS-CoV-2 y Minería de Datos.ipynb`: Libreta de Jupyter que nos ayudará a ir paso por paso por el taller.
- `SARSCoV2 y Minería de Datos.R`: Código fuente en lenguaje R que contiene las respuestas a la libreta de Jupyter.


Los primeros dos archivos fueron descargados del [National Center for Biotechnology Information (NCBI)](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Protein&VirusLineage_ss=Wuhan%20seafood%20market%20pneumonia%20virus,%20taxid:2697049&SLen_i=1273%20TO%201273&Completeness_s=complete&HostLineage_ss=Homo%20sapiens%20(human),%20taxid:9606) el día 7 de julio del 2021 y contienen datos de secuenciación de Estados Unidos de América.


#### Paqueterías necesarias

Para poder ejecutar correctamente el código, será necesaria la paquetería `seqinr`. Podemos instalarla de la isguiente manera:
```
install.packages("seqinr")
library(seqinr)
```

Base de datos adicional no utilizada en el taller: [GISAID Initiative](https://www.gisaid.org/)

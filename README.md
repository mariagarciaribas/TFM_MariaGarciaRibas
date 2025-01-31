#!/bin/sh

# Workflow para WGS


## 1. Control de calidad de las secuencias crudas mediante FastQC

fastqc *.fastqc.gz

multiqc .  # Para visulaizar todos los resultados de FastQC juntos


## 2. Limpieza y filtrado de secuencias con fastp

fastp \
    --thread 16 \  
    --verbose \  
    --in1 /PATH/TO/INPUT/SAMPLE_1.fq.gz \  # Archivo de entrada de lecturas forward
    --in2 /PATH/TO/INPUT/SAMPLE_2.fq.gz \  # Archivo de entrada de lecturas reverse
    --out1 /PATH/TO/OUTPUT/SAMPLE_fastp.1.fq.gz \  # Archivo de salida para lecturas forward después del procesamiento
    --out2 /PATH/TO/OUTPUT/SAMPLE_fastp.2.fq.gz \  # Archivo de salida para lecturas reverse después del procesamiento
    --unpaired1 /PATH/TO/OUTPUT/SAMPLE_fastp.unpaired.1.fq.gz \  # Archivo de salida para lecturas forward no pareadas
    --unpaired2 /PATH/TO/OUTPUT/SAMPLE_fastp.unpaired.2.fq.gz \  # Archivo de salida para lecturas reverse no pareadas
    --json /PATH/TO/OUTPUT/SAMPLE_report.json \  # Archivo JSON con el informe detallado de la ejecución
    --html /PATH/TO/OUTPUT/SAMPLE.html \  # Informe HTML visual de la calidad de las lecturas
    --report_title SAMPLE \  # Título del informe generado
    --detect_adapter_for_pe \  # Detecta y recorta adaptadores en lecturas emparejadas
    --correction \  # Corrige bases erróneas mediante consenso de secuencias
    --length_required 45 \  # Filtra lecturas con longitud menor a 45 bases
    --qualified_quality_phred 25 \  # Umbral de calidad mínima Phred para bases calificadas
    --average_qual 30 \  # Promedio mínimo de calidad Phred requerido en la lectura
    --cut_front \  # Recorta bases de baja calidad en el extremo 5'
    --cut_front_window_size 4 \  # Tamaño de la ventana para evaluar el recorte en el frente
    --cut_front_mean_quality 30 \  # Calidad media mínima para realizar el recorte en el frente
    --cut_tail \  # Recorta bases de baja calidad en el extremo 3'
    --cut_tail_window_size 4 \  # Tamaño de la ventana para evaluar el recorte en la cola
    --cut_tail_mean_quality 30 \  # Calidad media mínima para realizar el recorte en la cola
    --trim_front1 19 \  # Recorta los primeros 19 nucleótidos de las lecturas forward
    --trim_front2 19 \  # Recorta los primeros 19 nucleótidos de las lecturas reverse
    --trim_tail1 1 \  # Recorta 1 nucleótido del extremo 3' en lecturas forward
    --trim_tail2 1   # Recorta 1 nucleótido del extremo 3' en lecturas reverse

## 3. Detección de contaminación con Kraken y Bracken

kraken2 \  # Ejecuta Kraken2 para clasificar lecturas según su origen taxonómico
    --threads 8 \  # Usa 8 hilos para el análisis
    --db /PATH/TO/DB/kraken_std \  # Base de datos de referencia
    --fastq-input \  # Indica que las entradas son archivos FASTQ
    --unclassified-out unclassified_reads/SAMPLE.fq.gz \  # Salida de lecturas no clasificadas
    --classified-out classified_reads/SAMPLE.fq.gz \  # Salida de lecturas clasificadas
    --paired SAMPLE_1.fq.gz SAMPLE_2.fq.gz \  # Lecturas forward y reverse
    --report /PATH/TO/REPORTS/SAMPLE_kraken_report.txt \  # Informe de resultados de Kraken2
    --use-names  # Usa nombres de especies en el informe

bracken -d /PATH/TO/DB/kraken_std \  # Ejecuta Bracken para cuantificación de especies
    -i SAMPLE.kraken.report.txt \  # Archivo de entrada generado por Kraken2
    -o bracken.species.txt \  # Salida de abundancia de especies
    -l S  # Nivel taxonómico de especie

## 4. Ensamblaje de novo mediante SPAdes

spades.py --careful \  # Realiza un ensamblaje preciso
    -1 SAMPLE_clean_R1.fq \  # Lecturas forward limpias
    -2 SAMPLE_clean_R2.fq \  # Lecturas reverse limpias
    -o Spades_output  # Carpeta de salida

## 5. Evaluación de calidad del ensamblaje con QUAST

quast.py *.fasta -o quast  # Evalúa estadísticas de ensamblaje


## 6. Detección de contaminación en los genomas ensamblados con checkm

checkm lineage_wf \
    /PATH/TO/INPUT/scaffolds \  # Directorio con los ensamblajes en formato FASTA
    /PATH/TO/OUTPUT/checkm \  # Directorio de salida para los resultados de CheckM
    --extension fasta \  # Especifica que los archivos de entrada tienen extensión .fasta
    --file /PATH/TO/OUTPUT/checkm/Results.tsv \  # Archivo de salida con los resultados en formato tabular
    --tab_table \  # Genera una tabla con los resultados
    --threads 8   # Especifica el número de hilos para el procesamiento


## 7. Serotipificación 

### Serotipificación de los aislados de Salmonella con SeqSero2 y SISTR

SeqSero2_package.py \
    -p 10 \  # Número de hilos utilizados para el análisis
    -t 2 \  # Modo de análisis (2 para detección con datos de secuenciación en pares)
    -i /PATH/TO/INPUT/SAMPLE_1.fastq.gz /PATH/TO/INPUT/SAMPLE_2.fastq.gz  # Archivos de entrada de lecturas forward y reverse

sistr -f csv \
    -o /PATH/TO/OUTPUT/sistr/SAMPLE.csv \  # Archivo de salida en formato CSV con los resultados de serotipado
    /PATH/TO/INPUT/scaffolds/SAMPLE.fasta  # Archivo de entrada con el ensamblaje genómico

### Serotipificación de los aislados de Listeria monocytogenes con LisSero

lissero \
    --input /PATH/TO/INPUT/scaffolds/SAMPLE.fasta \  # Archivo de entrada con el ensamblaje genómico
    --output /PATH/TO/OUTPUT/lissero/SAMPLE.csv \  # Archivo de salida con los resultados del serotipado
    --threads 8  # Número de hilos utilizados para el análisis

## 8. Tipificación 

### Ariba para el análisis de MLST

ariba run \
    get_mlst/ref_db \  # Base de datos de referencia para MLST (Slamonella enterica y Listeria monocytogenes)
    /PATH/TO/INPUT/SAMPLE_1.fastq.gz \  # Archivo FASTQ con lecturas forward
    /PATH/TO/INPUT/SAMPLE_2.fastq.gz \  # Archivo FASTQ con lecturas reverse
    aribaMLST  # Directorio de salida con los resultados de MLST
#!/bin/bash

### chewBBACA para el análisis de cgMLST

#### Descargar esquema para Listeria monocytogenes
chewBBACA.py DownloadSchema \
    -sp 7 \  # Código de especie para Listeria monocytogenes
    -sc 1 \  # Nivel de confianza del esquema
    -o /home/usuari/cgMLST/Listeria_monocytogenes_Pasteur_cgMLST  # Directorio de salida

#### Descargar esquema para Salmonella enterica
chewBBACA.py DownloadSchema \
    -sp 8 \  # Código de especie para Salmonella enterica
    -sc 1 \  # Nivel de confianza del esquema
    -o /home/usuari/cgMLST/Salmonella_enterica_INNUENDO_cgMLST  # Directorio de salida

#### Llamada de alelos con chewBBACA
chewBBACA.py AlleleCall \
    -i /PATH/TO/INPUT/scaffolds \  # Directorio con ensamblajes en formato FASTA
    -g /home/usuari/cgMLST/Salmonella_enterica_INNUENDO_cgMLST or Listeria_monocytogenes_Pasteur_cgMLST \  # Esquema de cgMLST a utilizar
    -o /PATH/TO/OUTPUT/cgMLST \  # Directorio de salida
    --cpu 8  # Número de hilos utilizados para el análisis

### Construcción del árbol filogenético con Grapetreee

grapetree -p results_alleles.tsv -m MSTreeV2 >illumina.tree

## 9. Detección de elementos extracromosómicos (plásmidos)

### Ejecutar plasmidSPAdes para ensamblado de plásmidos

plasmidspades.py \
    --threads 8 \
    -o /PATH/TO/OUTPUT/plasmidspades \  # Directorio donde se almacenará el ensamblaje de plásmidos
    -1 /PATH/TO/INPUT/SAMPLE_1.fastq.gz \  # Lecturas forward de secuenciación
    -2 /PATH/TO/INPUT/SAMPLE_2.fastq.gz  # Lecturas reverse de secuenciación

### Ejecutar mob_recon para detección y caracterización de plásmidos

mob_recon \
    --infile /PATH/TO/INPUT/scaffolds/SAMPLE.fasta \  # Archivo de ensamblaje en formato FASTA
    --outdir /PATH/TO/OUTPUT/mob_recon  # Directorio de salida donde se almacenarán los resultados

### Ejecutar PlasmidFinder para identificar plásmidos basados en secuencias

plasmidfinder.py \
    -i /PATH/TO/INPUT/scaffolds/SAMPLE.fasta \  # Ensamblaje en formato FASTA
    -o /PATH/TO/OUTPUT/plasmidfinder/ \  # Directorio de salida para los resultados
    -p /PATH/TO/OUTPUT/plasmidfinder_db/ \  # Ruta de la base de datos de PlasmidFinder
    -x \  # Modo extendido de búsqueda
    -tmp /tmp/  # Carpeta temporal para el análisis

## 10. Detección de genes asociados a la viurlencia y a la resistencia a antimicrobianos (AMR)

### AMRFinderPlus para indentificar genes de AMR

amrfinder \
    -n /PATH/TO/INPUT/scaffolds/SAMPLE.fasta \  # Archivo de ensamblaje en formato FASTA
    --plus \  # Activa la detección extendida de genes de resistencia
    --threads 7 \  # Especifica el número de hilos para el procesamiento
    --name SAMPLE \  # Nombre de la muestra analizada
    -o /PATH/TO/OUTPUT/amrfinder/SAMPLE.tsv  # Archivo de salida con los resultados en formato TSV

### Abricate junto con la base de datos VFDB para detectar genes de virulencia 

abricate \
    --datadir /PATH/TO/abricate_db/ \  # Ruta de la base de datos de ABRicate
    --db VFDB \  # Base de datos seleccionada (VFDB para factores de virulencia)
    --threads 8 \  # Número de hilos utilizados en el análisis
    --minid 95 \  # Porcentaje mínimo de identidad para considerar una coincidencia
    --mincov 60 \  # Cobertura mínima requerida para la detección
    /PATH/TO/INPUT/scaffolds/SAMPLE.fasta \  # Ensamblaje genómico a analizar
    > /PATH/TO/OUTPUT/abricate/SAMPLE.tsv  # Archivo de salida con los resultados

## 11. Análisis de SNPs

### Identificación de especies para identificar el genoma de referencia  

kmerfinder.py \
    -i fastq.files \  # Directorio o archivo con las secuencias FASTQ de entrada
    -o /PATH/TO/OUTPUT/kmerfinder_results \  # Directorio de salida para los resultados
    -db /PATH/TO/DB/kmerfinder_db  # Ruta a la base de datos de KmerFinder

### Alineamiento filogenético con ParSNP

parsnp \
    -r REF_genomes/10403S.fna \  # Genoma de referencia para la alineación
    -d ST7 \  # Directorio con los genomas a analizar
    -o ST7/ST7_parSNP \  # Directorio de salida donde se almacenarán los resultados
    -p 7 \  # Número de hilos utilizados para acelerar el análisis
    -c \  # Habilita la alineación de genomas completos
    -x  # Excluye secuencias muy divergentes del análisis

### Cálculo de las distancia con snp-dists

snp-dists \
    test/good.aln \  # Archivo de alineamiento en formato FASTA generado por Parsnp u otro programa
    > distances.tab  # Archivo de salida con la matriz de distancias entre genomas


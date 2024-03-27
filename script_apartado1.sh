#!/bin/bash

mkdir -p Apartado1/output/fastqc


for sample in Apartado1/input/*.fastq
do
	echo "Ejecutando FastQC en $sample..."
	fastqc "$sample" -o Apartado1/output/fastqc
	echo "FastQC realizado para la muestra $sample"
done

mkdir -p Apartado1/output/index

# Indexar el genoma de referencia con Hisat2

for genome in Apartado1/input/*.fa
do
	genome_name=$(basename "$genome" .fa)
	index_name="Apartado1/output/index/${genome_name}_index"
	echo "Indexando el genoma $genome..."
	hisat2-build --seed 123 -p 2 "$genome" "$index_name"
done

echo "Indexación del genoma completada"

# Alineamiento con Hisat2 (en este caso vamos a usar unstranded por lo que no hace falta decir la orientacion, si no añadir --rna-strandness R)

mkdir -p Apartado1/output/hisat2

# Extracción del nombre del genoma indexado
genome_index_name=$(basename $(ls Apartado1/output/index/*.ht2 | head -1))
genome_index_name=${genome_index_name%.*.*}

## Hacemos un bucle a través del fastq de la lectura 1
for sample1 in Apartado1/input/*_1.fastq
do
	#Creamos variable para la lectura 2
	sample2="${sample1/_1.fastq/_2.fastq}"
	sample_name=$(basename "$sample1" _1.fastq)

	#Nombres de los outputs
	sam_output="Apartado1/output/hisat2/${sample_name}.sam"
	summary_output="Apartado1/output/hisat2/${sample_name}.hisat2.summary"

	echo "Alineando ${sample1} ..."
## Terminar ya que el sample name es solo para una muestra y aqui tenemos dos lecturas por muestra, entonces hay que cambiar eso (poner dos sample name si no) y modificarlo en el hisat2 para que cuadre

	hisat2 --new-summary --summary-file "$summary_output" \
	 --seed 123 --phred33 -p 2 -k 1 \
	-x Apartado1/output/index/"${genome_index_name}" \
	-1 "$sample1" -2 "$sample2" -S "$sam_output"
	echo "La muestra ${sample_name} ha sido alineada con exito"

done

echo "Alineamiento completado"

## Convertimos los SAM en BAM

for sam in Apartado1/output/hisat2/*.sam
do
	bam_name=$(basename "$sam" .sam)

	echo "Convirtiendo "$sam" a BAM ..."
	samtools view -bS "$sam" > "Apartado1/output/hisat2/${bam_name}.bam"

	echo "$sam convertido a BAM."
done


# Ordenamos los BAM
for bam in Apartado1/output/hisat2/*.bam
do
	bam_name_sorted=$(basename "$bam" .bam)
	echo "Ordenando $bam ..."
	samtools sort "$bam" -o "Apartado1/output/hisat2/${bam_name_sorted}.sorted.bam"

	echo "$bam ordenado"

done


## Indexamos los BAM ordenados

for sorted in Apartado1/output/hisat2/*.sorted.bam
do
	echo "Indexando $sorted"
	samtools index "$sorted"

done

echo "Todos los BAM indexados"



## Conteo de lecturas por gen mediante HTSEQ
echo "############################################################" \
	                                                            \
								    \

mkdir -p Apartado1/output/htseq
mkdir -p Apartado1/output/htseq/log

# Bucle para procesar cada archivo BAM ordenado
for bam in Apartado1/output/hisat2/*.sorted.bam
do
	name_bam=$(basename "$bam" .sorted.bam)

	echo "Procesando HTSeq-count para $name_bam"

	# Ejecutamos htseq-count (tampoco usamos --stranded=reverse, ya que hemos elegido unstranded)

	htseq-count --format=bam \
	--mode=intersection-nonempty \
	--minaqual=10 --type=exon --idattr=gene_id \
	--additional-attr=gene_name "$bam" \
	Apartado1/input/*.gtf > Apartado1/output/htseq/"$name_bam".htseq \
	2> Apartado1/output/htseq/log/"$name_bam"_htseq.log

	echo "HTSeq-count completado para $name_bam"
done


echo "#############################################################"

# Crear archivos BigWig

mkdir -p Apartado1/output/bigwig

for bam in Apartado1/output/hisat2/*.sorted.bam
do
	name_bam=$(basename "$bam" .sorted.bam)

	echo "Procesando BigWig para $name_bam"

	bamCoverage -b $bam -o Apartado1/output/bigwig/${name_bam}.bw --normalizeUsing BPM

# Comprobación del código de salida de bamCoverage

	if [ $? -eq 0 ]; then
		echo "$name_bam completado con éxito."
	else
		echo "Error al procesar $name_bam. bamCoverage falló."
		# Opcional: Salir del script si hay un error
		 exit 1
	fi

done

echo "Todos los archivos BigWig han sido creados con exito"



# Report MultiQC

echo "#############################################################"

echo "Generando MultiQC ... "
mkdir -p Apartado1/output/multiqc

multiqc -o Apartado1/output/multiqc Apartado1/output

if [ $? -eq 0 ]; then
	echo "Error en el procesamiento de MultiQC"
else
	echo "MultiQC completado"
fi



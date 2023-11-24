source /home/acampbe/software/miniconda3/bin/activate prokenv

annotated=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_SimulatedCC5CC8/refAnnotations/
mkdir -p $annotated
annotatedgffs=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_SimulatedCC5CC8/refGFFs/
mkdir -p $annotatedgffs

for file in /home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_SimulatedCC5CC8/refgenomes/* ; do

	basestring=$(basename $file)
	oldext=".fna"
	noext=""
	genomeid=${basestring/$oldext/$noext}
	outputdir=/home/acampbe/DFU/data/WGS_2020/cladebreaker/cladebreaker_SimulatedCC5CC8/annotations/$genomeid
	echo $genomeid
	prokka --outdir $outputdir --force --prefix $genomeid --genus Staphylococcus $file
	cp $outputdir/annot.faa $annotated$genomeid.faa
	cp $outputdir/annot.gff $annotatedgffs$genomeid.gff
done


#!/bin/bash
# Aggregate results from blast searches of crtOPQMN promoter in the assemblies 
# Into a single file

mkdir -p XanthinPromoterBlastResults
mv *.tab XanthinPromoterBlastResults/

for file in /home/acampbe/DFUStrainsWGS/StaphyloxanthinAlign/XanthinPromoterBlastResults/*.tab ; do

	base=$(basename $file)
        baseext="_cleaned.tab"
        blank=""
        extensionless=${base/$baseext/$blank}
	tab="       "
	prefix=$extensionless$tab
	sed -i -e "s/^/$prefix/" $file


done

cat /home/acampbe/DFUStrainsWGS/StaphyloxanthinAlign/XanthinPromoterBlastResults/*.tab > AllPromoterSequences.tab


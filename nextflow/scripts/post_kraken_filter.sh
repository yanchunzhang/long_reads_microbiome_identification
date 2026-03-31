prefix=$1
DB=$2
less $prefix.krakenuniq|awk '$1=="C"' > $prefix.krakenuniq.classified 
less $prefix.classified_by_krakenuniq |paste - - > $prefix.classified_by_krakenuniq.convert

krakenuniq-translate --db $DB --mpa-format $prefix.krakenuniq > $prefix.krakenuniq.translate
paste $prefix.classified_by_krakenuniq.convert $prefix.krakenuniq.translate $prefix.krakenuniq.classified \
	> $prefix.krakenuniq.info_collection && rm $prefix.classified_by_krakenuniq.convert $prefix.krakenuniq.translate $prefix.krakenuniq.classified
cat $prefix.krakenuniq.info_collection|awk '!/Eukaryota/ && /d__/ && !/synthetic/' \
	> $prefix.krakenuniq.info_collection.flt
less $prefix.krakenuniq.info_collection.flt|awk '{print ">"substr($1,2)"\n"$2}' > $prefix.krakenuniq.microbiome.fasta


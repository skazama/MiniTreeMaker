for i in `grep  "Now F" exe_16*.o* | cut -c 1-18`
do
	name=`echo $i | cut -c 5-15`
	echo "moving *${name}* to logfiles/"
	mv *${name}* logfiles/
done

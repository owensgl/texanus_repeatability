rep=$1

for selection in `seq 0 9`
do
	slim -d selection=0.$selection -d recom=1e-5 BC1_nonWF_underdominant.v0.2.slim | grep "^chr\|^01" > underdom/output.$rep.$selection.txt
	python3 ../cvtk/run_cvtk.py underdom/output.$rep.$selection.txt > underdom/conv.sel.$rep.$selection.txt
done

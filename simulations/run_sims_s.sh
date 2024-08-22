rep=$1

for selection in `seq 0 9`
do
	slim -d selection=0.$selection -d recom=1e-5 BC1_nonWF_bdm.v0.2.slim | grep "^chr\|^01" > output.$rep.$selection.txt
	python3 ../cvtk/run_cvtk.py output.$rep.$selection.txt > conv.sel.$rep.$selection.txt
done

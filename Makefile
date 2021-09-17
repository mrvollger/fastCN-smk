

fastCN/GC_control_gen:
	git clone https://github.com/KiddLab/fastCN.git \
		&& cd fastCN \
		&& sed -i -e 's/Use_strict 1/Use_strict 0/g' GC_control_gen.cc \
		&& g++ -o GC_control_gen GC_control_gen.cc \
		&& g++ -o SAM_GC_correction SAM_GC_correction.cc \
		&& cd .. \
		cargo install rustybam 


.PHONY: clean

clean:
	rm -rf fastCN


fastCN/GC_control_gen:
	which g++-11 \
		&& git clone https://github.com/KiddLab/fastCN.git \
		&& pushd fastCN \
		&& sed -i -e 's/Use_strict 1/Use_strict 0/g' GC_control_gen.cc \
		&& g++-11 -o GC_control_gen GC_control_gen.cc \
		&& g++-11 -o SAM_GC_correction SAM_GC_correction.cc \
		&& popd


.PHONY: clean

clean:
	rm -rf fastCN
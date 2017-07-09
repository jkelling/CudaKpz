all: bin effectiveExp slcg64calc combineTables selectValue smearToBin

bin: bin.cpp
	${CXX} --std=c++0x -o $@ $<
effectiveExp: effectiveExp.cpp
	${CXX} --std=c++0x -o $@ $< ../inc/libKMC.a
combineTables: combineTables.cpp
	${CXX} --std=c++0x -o $@ $< ../inc/libKMC.a
selectValue: selectValue.cpp
	${CXX} --std=c++0x -o $@ $< ../inc/libKMC.a
smearToBin: smearToBin.cpp
	${CXX} --std=c++11 -g -O0 -o $@ $< ../inc/libKMC.a

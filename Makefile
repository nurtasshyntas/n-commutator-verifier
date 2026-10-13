CXX ?= c++
CXXFLAGS ?= -O3 -std=c++17

.PHONY: all check verify verify-left verify-right clean
# Each rank-five run needs about 2.2 GiB; keep the large runs sequential.
.NOTPARALLEL:

all: direct_dense test_dense_literal

direct_dense: direct_dense.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

test_dense_literal: test_dense_literal.cpp direct_dense.cpp
	$(CXX) $(CXXFLAGS) test_dense_literal.cpp -o $@

check: test_dense_literal
	mkdir -p runs
	./test_dense_literal > runs/dense_literal_checks.json 2> runs/dense_literal_checks.log
	cat runs/dense_literal_checks.json

verify: verify-left verify-right

verify-left: direct_dense
	mkdir -p runs
	./direct_dense 5 Y left > runs/Y5_left.json 2> runs/Y5_left.log
	cat runs/Y5_left.json
	grep -F '"value":-3239792640,' runs/Y5_left.json > /dev/null

verify-right: direct_dense
	mkdir -p runs
	./direct_dense 5 Y right > runs/Y5_right.json 2> runs/Y5_right.log
	cat runs/Y5_right.json
	grep -F '"value":-3239792640,' runs/Y5_right.json > /dev/null

clean:
	rm -f direct_dense test_dense_literal

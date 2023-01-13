import pytest
import logging
import sys
print(sys.path)
logging.basicConfig()
logging.getLogger().setLevel(logging.DEBUG)

def test_matlab_output():
    import swash_pp.swash_mat2nc as snc
    path_run,run_file="src/swash_pp/tests/data/l31setup/","l31set03.sws"
    ds=snc.mat2nc(path_run=path_run,run_file=run_file)
    assert ds[0].Watlev.shape == (1,65,33601)

def test_table_output():
    import swash_pp.swash_mat2nc as snc
    path_run,run_file="src/swash_pp/tests/data/l31setup/","l31set03.sws"
    ds=snc.tab2nc(path_run=path_run,run_file=run_file)
    assert ds[0].Hsig.shape == (65,1)

if __name__ == '__main__':
    test_matlab_output()
    test_table_output()
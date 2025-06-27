#!/bin/bash
export PYTHONPATH=/home/zippy/vtkbool/build/lib/python3.13/site-packages/vtkbool:$PYTHONPATH
pushd testing

# pytest test_filter.py::test_triangle_strips -s

# pytest test_filter.py::test_equal_capt_pts -s
# pytest test_filter.py::test_equal_capt_pts_2 -s
# pytest test_filter.py::test_equal_capt_pts_3 -s

# gdb --args python -m pytest test_filter.py::test_simple -s

# pytest test_filter.py

pytest test_filter.py::test_simple -s

popd

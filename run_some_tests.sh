#!/bin/bash
export PYTHONPATH=/home/zippy/Desktop/vtkbool/build/lib/python3.14/site-packages/vtkbool:$PYTHONPATH
pushd testing

# pytest test_filter.py::test_simple -s
# gdb --args python -m pytest test_filter.py::test_simple -s

pytest test_filter.py

# ./generate_frieze.py

# pytest test_filter.py::test_simple -s
# pytest test_filter.py::test_simple_2 -s
# pytest test_filter.py::test_simple_3 -s
# pytest test_filter.py::test_simple_4 -s
# pytest test_filter.py::test_same -s
# pytest test_filter.py::test_intersecting_strips -s
# pytest test_filter.py::test_intersecting_strips_2 -s
# pytest test_filter.py::test_touch -s
# pytest test_filter.py::test_merger -s
# pytest test_filter.py::test_merger_2 -s
# pytest test_filter.py::test_quads -s
# pytest test_filter.py::test_triangle_strips -s
# pytest test_filter.py::test_special -s
# pytest test_filter.py::test_non_manifolds -s
# pytest test_filter.py::test_branched -s
# pytest test_filter.py::test_branched_2 -s
# pytest test_filter.py::test_branched_3 -s
# pytest test_filter.py::test_branched_4 -s
# pytest test_filter.py::test_branched_5 -s
# pytest test_filter.py::test_branched_6 -s
# pytest test_filter.py::test_self_intersecting_polys -s
# pytest test_filter.py::test_equal_capt_pts -s
# pytest test_filter.py::test_equal_capt_pts_2 -s
# pytest test_filter.py::test_equal_capt_pts_3 -s

popd

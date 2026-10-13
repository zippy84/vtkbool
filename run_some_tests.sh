#!/bin/bash
export PYTHONPATH=/home/zippy/vtkbool/build/lib/python3.14/site-packages/vtkbool:$PYTHONPATH
pushd testing

# pytest test_filter.py::test_simple -s
# gdb --args python -m pytest test_filter.py::test_simple -s

pytest test_filter.py

# ./generate_frieze.py

popd

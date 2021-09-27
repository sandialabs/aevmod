aevmod package
================
This package provides functionality for computing an atomic environment vector (AEV), as well as its Jacobian and Hessian. The AEV is a feature vector that is useful for representing the geometry of a molecule, or simply a set of atoms, in a manner that satisfies rotational and translational invariances. For the *i*-th atom in the system, its AEV is :math:`y_i=f_i(x)`, where *x* is a vector of Cartesian coordinates of all atoms in the set. The Jacobian and Hessian are with respect to the :math:`x` coordinates. Ther package uses `pybind11 <http://pybind11.readthedocs.io/en/stable/>`_ to expose our C++ AEV library to Python, as a python package ``aevmod``.  The package folder includes:

* this README file, which can be viewed in a web browser using `restview <https://pypi.python.org/pypi/restview>`_. You can use 

.. code-block:: shell 

    # simple invocation in present terminal
    restview README.rst &

    # or, invocation in a separate xterm terminal
    xterm -geometry 100x10-0-0 -e restview README.rst &

* a ``doc`` folder containing pdf documentation
* an ``examples`` folder containing an example application ``taev.py``
* an ``include`` folder containing requisite `sacado <https://trilinos.github.io/sacado.html>`_ include files for the C++ library
* a ``pyproject.toml`` file containing build-system dependencies and other configuration info
* a ``setup.py`` file for building and installing the package
* an ``src`` folder containing the C++ source code
* a ``tests`` folder containing code and data for running tests with `pytest <https://docs.pytest.org>`_ 

The ``aevmod`` package is compatible with Python 3, and the ``pybind11`` usage requires a C++ compiler that has at least C++11 functionality.

To try out ``aevmod``, under the package folder, do:

.. code-block:: shell

   # Build and install the package
   pip install .

   # Install pytest so we can run the tests
   pip install pytest

You can now run the tests to confirm everything works:

.. code-block:: shell

    pytest
    ====================================== test session starts =======================================
    platform darwin -- Python 3.8.10, pytest-6.2.4, py-1.10.0, pluggy-0.13.1
    rootdir: /Users/hnnajm/mls/pkg/aev, configfile: setup.cfg, testpaths: tests
    plugins: anyio-2.2.0
    collected 6 items

    tests/test_aev.py .                                                                        [ 16%]
    tests/test_ang_indsets.py .                                                                [ 33%]
    tests/test_hes.py .                                                                        [ 50%]
    tests/test_jac.py .                                                                        [ 66%]
    tests/test_rad_indsets.py .                                                                [ 83%]
    tests/test_structures.py .                                                                 [100%]

    ======================================= 6 passed in 0.75s ========================================

You can also run the example code:

.. code-block:: shell

    python ./examples/taev.py

If you make a code change you will need to re-build and install. You can
do this using pip:

.. code-block:: shell

    pip install --upgrade .




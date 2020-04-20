HOWTO
=====


Contribute to the documentation
-------------------------------

The documentation are located at `docs`.  If you want to help improving the
documentation you can do it simple as modifying/adding rst files inside the
`docs`. You need to have sphinx and the rtd them installed in your system 

.. code-block:: shell
   
   pip install -U sphinx
   pip install sphinx-rtd-theme

Once you modify the `docs` type

.. code-block:: shell
   
   cd docs 
   make html
   firefox index.html

You will need to reload the page once you run the `make`
again. 







.. include:: ../global.txt

.. _sec-install-arch:

Installing prerequisites and building PISM on Arch Linux
--------------------------------------------------------

PISM is registered as a `package <https://aur.archlinux.org/packages/pism/>`_
in the `Arch Linux User Repository <https://aur.archlinux.org>`_ (AUR). This package will build the latest
PISM release and its dependencies. PISM is linked to `Open MPI`_ and is
installed under ``/usr``.

Installation from the AUR requires an `AUR helper <https://wiki.archlinux.org/title/AUR_helpers>`_
such as `yay <https://aur.archlinux.org/packages/yay/>`_. If you do not already
have an AUR helper installed, install ``yay`` with the following commands:

.. code-block:: bash
   
   git clone https://aur.archlinux.org/yay.git
   cd git
   makepkg -si

You can then install PISM and its dependencies with the following command:

.. code-block:: bash
   
   yay -Sy pism

Once installed, the PISM binaries (e.g. ``pismr``, ``pismv``, various Python
tools) are available in the PATH and do not require further intervention to
work. It is recommended that the installation is manually verified with the
instructions in section :ref:`sec-install-quick-tests`.

.. note::

   The AUR package for PISM is maintained by Julien Seguinot; please e-mail him with
   questions about it. His email is: *first name two first letters + last name three first
   letters + at + posteo and then dot + eu*.

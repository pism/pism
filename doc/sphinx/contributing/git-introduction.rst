.. include:: ../global.txt

.. default-role:: literal

.. _sec-git-introduction:

Git introduction for PISM developers
====================================

.. _sec-git-configuration:

Recommended Git configuration
-----------------------------

Set name and e-mail address:

.. code-block:: bash

   git config --global user.name "John Doe"
   git config --global user.email johndoe@example.com

Do not push local branches nonexistent on upstream by default:

.. code-block:: bash

   git config --global push.default simple

Set an editor to use when writing commit messages. For example, run

.. code-block:: bash

   git config --global core.editor vi

to use `vi`.

.. _sec-new-feature-branch:

Working on a new branch
-----------------------

This section summarizes Git commands used in a typical development workflow. A good Git
GUI can save you from having to type these commands yourself but one should still know
them.

#. When starting work on a new feature make sure to start from the `dev` branch:

   .. code-block:: bash

      git checkout dev

   When working on a bug fix for the current (released) PISM version, start from `master`
   instead:

   .. code-block:: bash

      git checkout master

   See :ref:`sec-git-branches` for more.

#. Next, create and switch to a new branch:

   .. code-block::  bash

      git checkout -b <user-name>/<short-description>

   where `<user-name>` is your GitHub user name and `<short-description>` is a short (two
   or three words) description of the topic you intend to work on.

   For example:

   .. code-block:: bash

      git checkout -b ckhroulev/energy-balance

#. Work on the code, documentation, etc.

#. Inspect changes:

   .. code-block:: bash

      git status

#. Commit changes:

   - To commit all changes to files that are already in the repository:

     .. code-block:: bash

        git commit -a

   - To commit particular files

     .. code-block:: bash

        git commit file1 file2 ...

   - To add new files that are to be committed:

     .. code-block:: bash

        git add file1 file2 ...
        git commit

#. Push changes to GitHub:

   .. code-block:: bash

      git push -u origin ckhroulev/energy-balance

   Push changes to your own branch or other branches dedicated to a topic you may be
   sharing with your collaborators.

   .. note::

      Do *not* push to `dev` or `master` directly unless you know what you are doing.

#. If you started your branch from `dev` and need to use a feature that was added to `dev`
   since then you can run

   .. code-block:: bash

      git merge dev

   to "merge" the `dev` branch into your branch.

   .. note::

      Please do not merge `dev` into your branch unless you need to: doing it too often
      makes the development history (`git log`) more confusing and less useful.

.. _sec-better-commit-messages:

Writing better commit messages
------------------------------

Every commit should be accompanied by a meaningful message.

A commit message consists of a short one-line summary followed by a blank line and a
free-form body of the message.

The summary should be capitalized, use the imperative mood, and should not end in a period.

The body of a commit message should explain what changes were made and why (but not *how*).

If a commit addresses a GitHub issue, please include the issue title and number in the
body. Summarize the issue if its title is not descriptive enough.

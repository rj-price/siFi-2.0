"""Generated Qt resources for the GUI.

``ui_sifi2015.py``, ``ui_db_wizard.py`` and ``sifi_2015_rc.py`` are **generated**
from ``sifi2015.ui``, ``db_wizard.ui`` and ``sifi_2015.qrc`` beside them. Edit the
sources and regenerate; never hand-edit the outputs::

    cd src/sifi2/gui/resources
    pyrcc5 sifi_2015.qrc -o sifi_2015_rc.py
    pyuic5 --from-imports sifi2015.ui  -o ui_sifi2015.py
    pyuic5 --from-imports db_wizard.ui -o ui_db_wizard.py

``--from-imports`` matters. Without it pyuic5 ends its output with a bare
``import sifi_2015_rc``, which under Python 2's implicit relative imports found
the resource module next door but under Python 3 resolves against ``sys.path``.
Upstream had a 0-byte ``sifi_2015_rc.py`` at its repo root, so that import would
have succeeded, registered nothing, and rendered every icon blank
(``PLAN.md`` Phase 3, "resource-import hazard"). ``from . import sifi_2015_rc``
can only find the right module.

The ``.qrc`` lists five images, not upstream's twelve: three of the twelve
(``siFi21_logo_*.tif``) were never in the repository at all, four more
(``header0*.png``) are referenced by nothing, and the two TIFFs that *are* used
have been converted to PNG so Qt needs no TIFF image plugin.
"""

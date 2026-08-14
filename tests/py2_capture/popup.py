"""Stub for legacy/popup.py (PyQt4), so legacy/sifi_pipeline.py imports headlessly.

Placed ahead of legacy/ on sys.path by capture.py. Only the import needs to
succeed: capture.py never calls SifiPipeline.get_main_target, which is the sole
user of ListSelection.
"""


class ListSelection(object):
    def __init__(self, *args, **kwargs):
        raise NotImplementedError(
            "popup.ListSelection is stubbed out for fixture capture; "
            "get_main_target must not be called."
        )

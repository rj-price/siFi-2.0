"""Stub for legacy/show_plot.py (PyQt4 + matplotlib), so legacy/sifi_pipeline.py
imports headlessly.

Placed ahead of legacy/ on sys.path by capture.py. Only the import needs to
succeed: capture.py exercises the pure methods of SifiPipeline directly and never
calls run_pipeline, which is the sole user of DrawPlot.
"""


class DrawPlot(object):
    def __init__(self, *args, **kwargs):
        raise NotImplementedError(
            "show_plot.DrawPlot is stubbed out for fixture capture; "
            "run_pipeline must not be called."
        )

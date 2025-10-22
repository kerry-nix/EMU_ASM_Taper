"""
Init file for tools
"""

# from .automatic_bipolar_montage import automatic_bipolar_montage
# from .get_iEEG_data import get_iEEG_data
# from .gini import gini
# from .line_length import line_length
# from .pull_patient_localization import pull_patient_localization
# from .pull_sz_ends import pull_sz_ends
# from .pull_sz_starts import pull_sz_starts
# from .bandpower import bandpower
# from .movmean import movmean
# from .plot_iEEG_data import plot_iEEG_data
# from .clean_labels import clean_labels
# from .find_non_ieeg import find_non_ieeg
import pkgutil
import importlib
import tools   # assuming X is the name of the package

# iterate over all module names in the X package
for loader, name, is_pkg in pkgutil.walk_packages(tools.__path__):
    # import the module and add its functions to the package namespace
    module = importlib.import_module(f"{tools.__name__}.{name}")
    for key, value in module.__dict__.items():
        if callable(value):
            setattr(tools, key, value)

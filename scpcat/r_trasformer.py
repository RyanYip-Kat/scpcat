import os
from functools import wraps
from typing import Optional



class MissingR(Exception):
    def __init__(self):
        super().__init__('Missing R setup in current system')


class RRuntimeException(Exception):
    def __init__(self, msg):
        super().__init__('R Runtime Exception: {}'.format(msg))

class MissingFunctionException(Exception):
    def __init__(self):
        super().__init__('Selected function not available')

def ensure_R_setup():
    from rpy2 import situation
    try:
        if not situation.get_r_home() or not situation.r_version_from_subprocess():
            raise MissingR()

    except MissingR as e:
        raise e


def _ensure_path_exists(path: str) -> None:
    expanded_path = os.path.expanduser(path)

    if not os.path.exists(expanded_path):
        os.makedirs(expanded_path)


def with_r_setup(f):
    @wraps(f)
    def wrapper(*args, **kwargs):
        ensure_R_setup()

        from rpy2.rinterface_lib.embedded import RRuntimeError
        from rpy2 import robjects

        return f(*args, **kwargs, robjects=robjects, r_runtime_error=RRuntimeError)

    return wrapper


@with_r_setup
def Seurat2AnnData(*,
        filename:str,
        name:str,
        outdir:str,
        robjects,
        r_runtime_error: Exception) -> None:
    this_file_dir = os.path.dirname(os.path.realpath(__file__))
    robjects.r.source(os.path.join(this_file_dir, 'R/Seurat2Anndata.R'))
    available_names = list(robjects.globalenv.keys())
    main_function = 'Seurat2AnnData'
    
    if main_function in available_names:
        function_name = main_function
    else:
        raise MissingFunctionException()

    transformer = robjects.r[function_name]
    try:
        transformer(filename=filename,name=name,outDir=outdir)
    except r_runtime_error as e:
        raise RRuntimeException(e)

        


from pathlib import Path
from typing import Union
import logging

logger = logging.getLogger("metagraph_workflow")


def create_transcript_path_list(path: Union[Path, str], transcript_path: Union[Path, str], suffix=''):
    paths = [str(p.absolute()) for p in Path(path).glob(f'*{suffix}')]

    with open(transcript_path, 'w') as f:
        f.write('\n'.join(paths))

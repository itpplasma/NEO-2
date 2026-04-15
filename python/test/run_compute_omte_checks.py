"""Run compute_omte regression checks without requiring pytest."""

from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))


def main():
    module_path = ROOT / 'test' / 'test_compute_omte.py'
    spec = spec_from_file_location('test_compute_omte', module_path)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    tests = [
        getattr(module, name)
        for name in sorted(dir(module))
        if name.startswith('test_') and callable(getattr(module, name))
    ]
    for test in tests:
        test()
    print(f'Ran {len(tests)} compute_omte checks')


if __name__ == '__main__':
    main()

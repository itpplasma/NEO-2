from pathlib import Path

import numpy as np


def create_test_mars_dir(base_dir):
    base_dir = Path(base_dir)
    base_dir.mkdir(parents=True, exist_ok=True)

    sqrtspol = np.linspace(0.0, 1.0, 11)
    sqrtstor = np.sqrt(sqrtspol)
    n_e = 1.0e19 * (1.0 + sqrtspol)
    t_e = 1.0e3 * (1.0 + 0.5 * sqrtspol)
    t_i = 1.5e3 * (1.0 + 0.25 * sqrtspol)
    vrot = 1.0e4 * sqrtspol

    _write_profile(base_dir / 'PROFDEN.IN', sqrtspol, n_e)
    _write_profile(base_dir / 'PROFTE.IN', sqrtspol, t_e)
    _write_profile(base_dir / 'PROFTI.IN', sqrtspol, t_i)
    _write_profile(base_dir / 'PROFROT.IN', sqrtspol, vrot)
    np.savetxt(base_dir / 'sqrtstor.dat', np.column_stack([sqrtspol, sqrtstor]))

    (base_dir / 'RUN.IN').write_text(
        '&KINETIC\n'
        ' ESPECIES_Z = -1, 1\n'
        ' ESPECIES_M = 0.000544617, 2.0\n'
        '/\n',
        encoding='utf-8',
    )

    return {
        'mars_dir': base_dir,
        'sqrtspol': sqrtspol,
        'sqrtstor': sqrtstor,
        'n_e': n_e,
        'n_i': n_e,
        't_e': t_e,
        't_i': t_i,
        'vrot': vrot,
    }


def _write_profile(path, x, y):
    header = f'{len(x)} 1'
    np.savetxt(path, np.column_stack([x, y]), header=header, comments='')

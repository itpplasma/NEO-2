#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np

def write_profiles_to_dat_files(profiles_mars: dict, output_dir: str):
    os.makedirs(output_dir, exist_ok=True)
    for profile_name, profile_data in profiles_mars.items():
        filename = f"{output_dir}/{profile_name}.dat"
        np.savetxt(filename, profile_data)
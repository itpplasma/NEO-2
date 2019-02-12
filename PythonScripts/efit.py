#!/usr/bin/env python3

def append_array_to_string_in_five_columns(lines: str, array):
  nr_elements = len(array)

  formatstring_five_elements = ' {0[0]: 1.8e} {0[1]: 1.8e} {0[2]: 1.8e} {0[3]: 1.8e} {0[4]: 1.8e}\n'

  for i in range(int(nr_elements/5)):
    lines += formatstring_five_elements.format(array[i*5:i*5+5])

  if (nr_elements % 5 == 1):
    lines += ' {0[0]: 1.8e}'.format(array[-1:])
  elif (nr_elements % 5 == 2):
    lines += ' {0[0]: 1.8e} {0[1]: 1.8e}'.format(array[-2:])
  elif (nr_elements % 5 == 3):
    lines += ' {0[0]: 1.8e} {0[1]: 1.8e} {0[2]: 1.8e}'.format(array[-3:])
  elif (nr_elements % 5 == 4):
    lines += ' {0[0]: 1.8e} {0[1]: 1.8e} {0[2]: 1.8e} {0[3]: 1.8e}'.format(array[-4:])

  return lines

"""Extract and store information from an efit file.

This class is designed to extract and store information from an efit
file (textformat, but not really human readable).

The constructor/init-method expects a file name, from these the
information are read.
"""
class efit_type:
  characters_per_number=16
  numbers_per_line=5

  def __init__(self, filename: str):
    import numpy as np

    # Read in the file.
    with open(filename) as f:
        lines = f.readlines()

    # Split the first line, which contains some usefull (e.g. nw, nh) and
    # some not so usefull information (e.g. unknownvaluestr and signtodiscard).
    self.unknownvaluestr, self.datestr, self.signtodiscard, self.shootnr, self.timestr, self.idum,self.nw,self.nh = lines[0].split()

    # Variables are so far all strings, but change to integer for some of them.
    self.shootnr = int(self.shootnr)
    self.idum = int(self.idum)# Some
    self.nw = int(self.nw)# Number of horizontal R grid points
    self.nh = int(self.nh)# Number of vertical z grid points

    # First few (four) lines contain specific values.
    self.rdim, self.zdim, self.rcentr, self.rleft, self.zmid = np.array(lines[1].split(),dtype=float)
    self.rmaxis, self.zmaxis, self.simag, self.sibry, self.bcentr = np.array(lines[2].split(),dtype=float)
    self.current, self.simag2, self.dummy1, self.rmaxis, self.dummy2 = np.array(lines[3].split(),dtype=float)
    self.zmaxis, self.dummy3, self.sibry2, self.dummy4, self.dummy5 = np.array(lines[4].split(),dtype=float)
    #~ self.rmaxis.unit = "m"
    #~ self.zmaxis.unit = "m"
    #~ self.simag.unit = "Weber/rad"
    #~ self.sibry.unit = "Weber/rad"

    # Create array of psi points.
    self.psi = np.linspace(0,self.sibry-self.simag,self.nw)

    # Now read in six array, the first and the last one of size nw, the fourth of size nw*nh.
    offset = 5
    self.fpol = self.local_split(''.join(lines[offset:offset+int(self.nw/self.numbers_per_line)]))

    offset = offset+int(self.nw/self.numbers_per_line)
    self.pres = self.local_split(''.join(lines[offset:offset+int(self.nw/self.numbers_per_line)]))
    #~ self.pres.unit = "N/m^2"

    offset = offset+int(self.nw/self.numbers_per_line)
    self.ffprim = self.local_split(''.join(lines[offset:offset+int(self.nw/self.numbers_per_line)]))
    #~ self.ffprim.unit = "(mT)^2/(Weber/rad)"

    offset = offset+int(self.nw/self.numbers_per_line)
    self.pprime = self.local_split(''.join(lines[offset:offset+int(self.nw/self.numbers_per_line)]))
    #~ self.pprime.unit = "(N/m^2)/(Weber/rad)"

    offset = offset+int(self.nw/self.numbers_per_line)
    self.psi2d = self.local_split(''.join(lines[offset:offset+int(self.nw*self.nh/self.numbers_per_line)])).reshape([self.nh,self.nw])
    #~ self.psi2d.unit = "Weber/rad"

    offset = offset+int(self.nw*self.nh/self.numbers_per_line)
    self.q = self.local_split(''.join(lines[offset:offset+int(self.nw/self.numbers_per_line)]))

    # Next are two sizes in the file.
    offset = offset+int(self.nw/self.numbers_per_line)
    self.number_boundary_points, self.number_limiter_points = lines[offset].split()
    self.number_boundary_points = int(self.number_boundary_points)
    self.number_limiter_points = int(self.number_limiter_points)

    # It follows an array of (r,z) values of the boundary points.
    offset = offset+1
    self.unknown_array1 = self.local_split(''.join(lines[offset:offset+int(np.ceil(2*self.number_boundary_points/self.numbers_per_line))])).reshape([self.number_boundary_points,2])

    # The following code does not work, as the format of the file is not
    # as expected.
    # For the default file and another one there are 4*5 additional
    # values (0) at the end of the array.
    offset = offset+int(2*self.number_boundary_points/self.numbers_per_line)
    # \attention The '+10' is due to the actual layout of (some?) files. It is not sure wether this is correct or not.
    self.unknown_array2 = self.local_split(''.join(lines[offset:offset+int(np.ceil(2*(self.number_limiter_points+10)/self.numbers_per_line))])).reshape([self.number_limiter_points,2])

    # Rest of the the file has so far unknown content. There seems to
    # follow a line with three values, another array and
    # a last value?
    # If it is known what this is, it is read, until then we just ignore
    # it.

  def write(self, filename: str):
    import numpy as np

    lines = ''
    formatstring_five_elements = ' {: 1.8e} {: 1.8e} {: 1.8e} {: 1.8e} {: 1.8e}'

    # Split the first line, which contains some usefull (e.g. nw, nh) and
    # some not so usefull information (e.g. unknownvaluestr and signtodiscard).
    lines += '   ' + self.unknownvaluestr + '   ' + self.datestr + '      ' + self.signtodiscard + ' {}  '.format(self.shootnr) + self.timestr + '           {: >3d}'.format(self.idum) + ' {: >3d}'.format(self.nw) + ' {: >3d}'.format(self.nh) + '\n'

    # First few (four) lines contain specific values.
    lines += formatstring_five_elements.format(self.rdim, self.zdim, self.rcentr, self.rleft, self.zmid) + '\n'
    lines += formatstring_five_elements.format(self.rmaxis, self.zmaxis, self.simag, self.sibry, self.bcentr) + '\n'
    lines += formatstring_five_elements.format(self.current, self.simag2, self.dummy1, self.rmaxis, self.dummy2) + '\n'
    lines += formatstring_five_elements.format(self.zmaxis, self.dummy3, self.sibry2, self.dummy4, self.dummy5) + '\n'

    lines = append_array_to_string_in_five_columns(lines, self.fpol)
    lines = append_array_to_string_in_five_columns(lines, self.pres)
    lines = append_array_to_string_in_five_columns(lines, self.ffprim)
    lines = append_array_to_string_in_five_columns(lines, self.pprime)

    lines = append_array_to_string_in_five_columns(lines, np.ravel(self.psi2d))

    lines = append_array_to_string_in_five_columns(lines, self.q)

    lines += ' {: >4d} {: >4d}\n'.format(self.number_boundary_points, self.number_limiter_points)

    lines = append_array_to_string_in_five_columns(lines, np.ravel(self.unknown_array1))

    lines = append_array_to_string_in_five_columns(lines, np.ravel(self.unknown_array2))

    lines += '\n'

    with open(filename, 'w') as f:
      f.writelines(lines)

  def local_split(self, a):
    import numpy as np

    aflat = a.replace('\n','')
    alist = [aflat[self.characters_per_number*k:self.characters_per_number*(k+1)] for k in range(int(len(aflat)/self.characters_per_number))]
    return np.array(alist, dtype=float)

"""Function to compute isolines for drawing.

Algorythm taken from http://paulbourke.net/papers/conrec/ (c-version to
be specific).
"""
def conrec(data, nr_horizontal_points, nr_vertical_points, contourlevels, xcoordinates, ycoordinates):
  im = [0,1,1,0]
  jm = [0,0,1,1]
  case_table = [
     [ [0,0,8],[0,2,5],[7,6,9] ],
     [ [0,3,4],[1,3,1],[4,3,0] ],
     [ [9,6,7],[5,2,0],[8,0,0] ]
   ]
  h = [0.0, 0.0, 0.0, 0.0, 0.0]
  sh = [0, 0, 0, 0, 0]
  xh = [0.0, 0.0, 0.0, 0.0, 0.0]
  yh = [0.0, 0.0, 0.0, 0.0, 0.0]

  x1 = 0.0
  x2 = 0.0
  y1 = 0.0
  y2 = 0.0

  # Local functions.
  def Line_between_vertices_1_and_2():
    x1 = xh[m1]
    y1 = yh[m1]
    x2 = xh[m2]
    y2 = yh[m2]
    return [x1, y1, x2, y2]

  def Line_between_vertices_2_and_3():
    x1 = xh[m2]
    y1 = yh[m2]
    x2 = xh[m3]
    y2 = yh[m3]
    return [x1, y1, x2, y2]

  def Line_between_vertices_3_and_1():
    x1 = xh[m3]
    y1 = yh[m3]
    x2 = xh[m1]
    y2 = yh[m1]
    return [x1, y1, x2, y2]

  def Line_between_vertex_1_and_side_2_3():
    x1 = xh[m1]
    y1 = yh[m1]
    x2 = xsect(m2,m3)
    y2 = ysect(m2,m3)
    return [x1, y1, x2, y2]

  def Line_between_vertex_2_and_side_3_1():
    x1 = xh[m2]
    y1 = yh[m2]
    x2 = xsect(m3,m1)
    y2 = ysect(m3,m1)
    return [x1, y1, x2, y2]

  def Line_between_vertex_3_and_side_1_2():
    x1 = xh[m3]
    y1 = yh[m3]
    x2 = xsect(m1,m2)
    y2 = ysect(m1,m2)
    return [x1, y1, x2, y2]

  def Line_between_sides_1_2_and_2_3():
    x1 = xsect(m1,m2)
    y1 = ysect(m1,m2)
    x2 = xsect(m2,m3)
    y2 = ysect(m2,m3)
    return [x1, y1, x2, y2]

  def Line_between_sides_2_3_and_3_1():
    x1 = xsect(m2,m3)
    y1 = ysect(m2,m3)
    x2 = xsect(m3,m1)
    y2 = ysect(m3,m1)
    return [x1, y1, x2, y2]

  def Line_between_sides_3_1_and_1_2():
    x1 = xsect(m3,m1)
    y1 = ysect(m3,m1)
    x2 = xsect(m1,m2)
    y2 = ysect(m1,m2)
    return [x1, y1, x2, y2]

  def xsect(p1, p2):
    return ((h[p2]*xh[p1]-h[p1]*xh[p2])/(h[p2]-h[p1]))

  def ysect(p1, p2):
    return ((h[p2]*yh[p1]-h[p1]*yh[p2])/(h[p2]-h[p1]))

  options = {
      1 : Line_between_vertices_1_and_2,
      2 : Line_between_vertices_2_and_3,
      3 : Line_between_vertices_3_and_1,
      4 : Line_between_vertex_1_and_side_2_3,
      5 : Line_between_vertex_2_and_side_3_1,
      6 : Line_between_vertex_3_and_side_1_2,
      7 : Line_between_sides_1_2_and_2_3,
      8 : Line_between_sides_2_3_and_3_1,
      9 : Line_between_sides_3_1_and_1_2}

  contourlines = []
  for k in range(len(contourlevels)):
    contourlines.append([])

  for j in range(nr_vertical_points-2, 0-1, -1):
    for i in range(0, nr_horizontal_points-1, 1):
      dmin = min(data[j][i], data[j+1][i], data[j][i+1], data[j+1][i+1])
      dmax = max(data[j][i], data[j+1][i], data[j][i+1], data[j+1][i+1])
      if (dmax < contourlevels[0] or dmin > contourlevels[-1]):
        continue
      for k in range(len(contourlevels)):
        if (dmax < contourlevels[k] or dmin > contourlevels[k]):
          continue
        for m in range(5-1, 0-1, -1):
          if m > 0:
            h[m]  = data[j+jm[m-1]][i+im[m-1]]-contourlevels[k]
            xh[m] = float(xcoordinates[i+im[m-1]])
            yh[m] = float(ycoordinates[j+jm[m-1]])
          else:
            h[0]  = 0.25 * float(h[1]+h[2]+h[3]+h[4])
            xh[0] = 0.50 * float(xcoordinates[i]+xcoordinates[i+1])
            yh[0] = 0.50 * float(ycoordinates[j]+ycoordinates[j+1])
            if (h[m] > 0.0):
              sh[m] = 1
            elif (h[m] < 0.0):
              sh[m] = -1
            else:
              sh[m] = 0

        # Note: at this stage the relative heights of the corners and the
        # centre are in the h array, and the corresponding coordinates are
        # in the xh and yh arrays. The centre of the box is indexed by 0
        # and the 4 corners by 1 to 4 as shown below.
        # Each triangle is then indexed by the parameter m, and the 3
        # vertices of each triangle are indexed by parameters m1,m2,and m3.
        # It is assumed that the centre of the box is always vertex 2
        # though this isimportant only when all 3 vertices lie exactly on
        # the same contour level, in which case only the side of the box
        # is drawn.
        #   vertex 4 +-------------------+ vertex 3
        #            | \               / |
        #            |   \    m-3    /   |
        #            |     \       /     |
        #            |       \   /       |
        #            |  m=2    X   m=2   |       the centre is vertex 0
        #            |       /   \       |
        #            |     /       \     |
        #            |   /    m=1    \   |
        #            | /               \ |
        #   vertex 1 +-------------------+ vertex 2

        # Scan each triangle in the box
        for m in range(1, 5, +1):
          m1 = m
          m2 = 0
          if (m != 4):
            m3 = m + 1
          else:
            m3 = 1
          case_value = case_table[sh[m1]+1][sh[m2]+1][sh[m3]+1]
          if (case_value == 0):
            continue
          #~ https://bytebaker.com/2008/11/03/switch-case-statement-in-python/
          try:
            contourlines[k].append(options[case_value]())
          except KeyError:
            print('Ignoring problem in selecting the right case.')

  return contourlines

def compute_rgrid(efitvar):
  # \todo Check for correct computation.
  rgrid = []
  for i in range(efitvar.nw):
    rgrid.append(i)
  return rgrid

def compute_zgrid(efitvar):
  # \todo Check for correct computation.
  zgrid = []
  for i in range(efitvar.nh):
    zgrid.append(i)
  return zgrid

def compute_contourlevels(efitvar, number_of_fluxsurfaces):
  # \todo Check for correct computation.
  contourgrid = []
  for i in range(number_of_fluxsurfaces):
    contourgrid.append(efitvar.simag + (efitvar.sibry - efitvar.simag)/(number_of_fluxsurfaces-1)*i)
  return contourgrid

"""Check if a given path is closed.

It is assumed to be closed, if the distance between the last and the
first point is not much greater than the distance between the first and
the second point. 'Not much greater' is defined by an internal consant.
"""
def closed(p):
  acceptance_factor = 1.5
  #~ v = p.vertices
  #~ x = v[:,0]
  #~ y = v[:,1]
  x = []
  y = []
  for v in p.iter_segments():
    x.append(v[0][0])
    y.append(v[0][1])
  # Maybe change to comparison with average over whole range.
  if (x[0] - x[-1])**2 + (y[0] - y[-1])**2 < acceptance_factor*((x[0] - x[1])**2 + (y[0] - y[1])**2):
    return True
  else:
    return False

"""Takes function with discretization and recomputes the values for another discretization.

This function takes the values of a function at given positions, and will
calculate the values at different, given positions.
It is assumed that the new grid is inside the old one, i.e. only
interpolation of the data is required, not extrapolation.
"""
def recompute_from_to(var_on_ingrid, ingrid, outgrid):
  import numpy as np
  from scipy.interpolate import interp1d

  f = interp1d(ingrid, var_on_ingrid)
  var_on_outgrid = f(outgrid)

  return var_on_outgrid

"""Make sure outputs x and y have given length.

This function takes two lists x,y with the same number of elements, and
output two arrays with 'length' number of elements, assuming they are
periodic.
This is done by using function 'recompute_from_to' with appropriate
input and output grids.
"""
def make_length(x, y, length):
  import numpy as np

  if (len(x) != len(y)):
    raise Exception
  else:
    ingrid = np.array(range(len(x)+1), dtype=float)/float(len(x))
    outgrid = np.array(range(length), dtype=float)/float(length)
    var_on_outgrid_x = recompute_from_to(x + [x[0]], ingrid, outgrid)
    var_on_outgrid_y = recompute_from_to(y + [y[0]], ingrid, outgrid)

    return [var_on_outgrid_x, var_on_outgrid_y]

def make_fouriertrafo(x):
  import numpy.fft as fft

  return fft.rfft(x) #[1:]

"""Expects filename of a efit file, and will print a boozer file with given name.

This function has two names as input arguments. The first one is the
name of an efit file to be read and converted to a boozer file. The
result is writen in a file with the second name.
"""
def efit_to_boozer(filename, outputfile):
  import matplotlib.pyplot as plt
  import numpy as np

  import boozer

  number_of_fluxsurfaces=30
  efitvar = efit_type(filename)

  rgrid = compute_rgrid(efitvar)
  zgrid = compute_zgrid(efitvar)

  #~ contourlevels = compute_contourlevels(efitvar, number_of_fluxsurfaces)
  #~ contourlines = conrec(efitvar.psi2d, efitvar.nw, efitvar.nh, contourlevels, rgrid, zgrid)
  contourlines = []
  cn = plt.contour(efitvar.psi2d,number_of_fluxsurfaces)
  contourlevels = cn.levels

  # Determine the number of closed surfaces in advance, as this is equal
  # to the parameter nsurf, required for the header.
  j = 0
  for k in range(len(cn.collections)):
    for p in cn.collections[k].get_paths():
      if (closed(p)):
        j = j+1

  # \todo Use correct values.
  m0b =20 #
  n0b = 0 # Data only for single poloidal cut, thus only zero mode can be calculated.
  nsurf = j
  nfp = 1 #
  psi_tor_a = -3.0 #
  aminor = 0.5 #
  Rmajor = efitvar.rmaxis
  vmnb = np.zeros(m0b+1) #
  bmnb = np.zeros(m0b+1) #
  radial_restriction_parameter = 1.01
  s = (cn.levels - radial_restriction_parameter*min(cn.levels))/(radial_restriction_parameter*(max(cn.levels) - min(cn.levels))) #
  iota_ = 1.0/efitvar.q
  # Scale efitvar.psi to be in the range [0,1] as s is scaled that way.
  iotavec = recompute_from_to(iota_, efitvar.psi/max(efitvar.psi), s)
  bsubvB = 0.0 #
  bsubuB = 0.0 #
  pprime = 0.0 #
  vp = 0.0 #
  nb = np.zeros(m0b+1, dtype=int)
  mb = range(m0b+1)
  shot = efitvar.shootnr

  boozer.write_boozer_head(outputfile, '01', shot, m0b, n0b, nsurf, nfp, psi_tor_a, aminor, Rmajor)

  for (k, psi, iota) in zip(range(len(cn.collections)-1,0-1, -1), s, iotavec):
    for p in cn.collections[k].get_paths():

      if (closed(p)):
        x = []
        z = []
        for v in p.iter_segments():
          x.append(v[0][0])
          z.append(v[0][1])

        # Final length chosen such, that there will be m0b modes, not counting the zero mode (constant) and the first order (sine(x)).
        [xl, zl] = make_length(x, z, 2*(m0b + 0))
        # Make the fouriertrafo.
        xf = make_fouriertrafo(xl)
        zf = make_fouriertrafo(zl)

        boozer.append_boozer_block_head(outputfile, psi, iota, bsubvB, bsubuB, pprime, vp, nfp)
        boozer.append_boozer_block(outputfile, mb, nb, xf, zf, vmnb, bmnb, nfp)

  # As test plot the contour lines.
  #~ plt.figure()
  #~ colormap = plt.get_cmap('cool', 30)
  #~ for k in range(number_of_fluxsurfaces):
    #~ for l in range(len(contourlines[k])):
      #~ plt.plot([contourlines[k][l][0], contourlines[k][l][2]], [contourlines[k][l][1], contourlines[k][l][3]], color=colormap(float(k)/float(number_of_fluxsurfaces-1)))
  #~ plt.savefig('test')

########################################################################
if __name__ == "__main__":
  import matplotlib.pyplot as plt
  import sys

  # Check input arguments for input and output filename.
  if (len(sys.argv) >= 2):
    filename = sys.argv[1]
  else:
    filename = '/proj/plasma/Neo2/ASDEX-U/31114/g031114.05000'

  if (len(sys.argv) >= 3):
    outputfile = sys.argv[2]
  else:
    outputfile = 'boozer.bc'

  efit_to_boozer(filename, outputfile)

  #~ efitvar = efit_type(filename)

  #~ plt.figure()
  #~ plt.plot(efitvar.psi,efitvar.fpol)
  #~ plt.xlabel('psi_pol')
  #~ plt.ylabel('F_pol')
  #~ plt.savefig('F_pol')

  #~ plt.figure()
  #~ plt.plot(efitvar.psi,efitvar.pres)
  #~ plt.xlabel('psi_pol')
  #~ plt.ylabel('pressure')
  #~ plt.savefig('pressure')

  #~ plt.figure()
  #~ plt.plot(efitvar.psi,efitvar.ffprim)
  #~ plt.xlabel('psi_pol')
  #~ plt.ylabel('F*dF/dpsi_pol')
  #~ plt.savefig('FdF_dpsi_pol')

  #~ plt.figure()
  #~ plt.plot(efitvar.psi,efitvar.pprime)
  #~ plt.xlabel('psi_pol')
  #~ plt.ylabel('dp/dpsi_pol')
  #~ plt.savefig('dp_dpsi_pol')

  #~ plt.figure()
  #~ plt.contour(efitvar.psi2d,30)
  #~ plt.colorbar()
  #~ plt.axis('equal')
  #~ plt.savefig('equal')

  #~ plt.figure()
  #~ plt.plot(efitvar.psi,efitvar.q)
  #~ plt.ylim(0,5)
  #~ plt.xlabel('psi_pol')
  #~ plt.ylabel('safety factor q')
  #~ plt.savefig('safety_factor')
